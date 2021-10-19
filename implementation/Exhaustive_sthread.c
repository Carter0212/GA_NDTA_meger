#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "node.h"

//------------structure delcare------------
struct max_ue_record ue_rec[ue_num];
struct max_bs_record bs_rec;
struct ue_buffer ue_buf[4][ue_num];
struct ue_parameters ue_ind[ue_num];
struct cbsd_parameters cbsd_ind[cbsd_num];
struct building_parameters building_ind[cbsd_num][2];
struct mbs_parameters mbs_ind;

//------------output result------------
void write_excel(){ // write the max average power times result to excel 
	int i;
	
	FILE *fpw=NULL;
	fpw=fopen("Exhaustive.csv","w");
	fprintf(fpw,"the max average throughput times\n\n");
	
	fprintf(fpw,"BS TX power\n");
	for(i=0;i<4;i++) fprintf(fpw,"MBS area %d,",i+1);
	for(i=0;i<cbsd_num;i++) fprintf(fpw,"CBSD %d,",i+1);
	fprintf(fpw,"\n");
	for(i=0;i<4;i++) fprintf(fpw,"%.3f,",bs_rec.a_p[i]);
	for(i=0;i<cbsd_num;i++) fprintf(fpw,"%.3f,",bs_rec.c_p[i]);
	fprintf(fpw,"\n");

	fprintf(fpw,"max average power\n");
	fprintf(fpw,"%f\n\n",bs_rec.avg_t);	
	
	fprintf(fpw,"x,y,area,be blocked (1 : NLOS / 0 : LOS),be blocked number (>=1 : CBSD(x) / 0 : no),connect bs number (>=1 : CBSD / 0 : MBS),throughput (unit : M bps)\n");
	for(i=0;i<ue_num;i++) fprintf(fpw,"%d,%d,%d,%d,%d,%d,%f\n",ue_rec[i].x,ue_rec[i].y,ue_rec[i].area,ue_rec[i].block,ue_rec[i].block_num,ue_rec[i].connect_bs_num,ue_rec[i].throughput);
	
	fclose(fpw);
}

//------------choose higher bs connect------------
int choose_better_bs(int area, int i, int ue_c){
	float best_bs; // record the best bs number
	int best_c; // for count
	
	// choose the better bs to connect
	best_bs=ue_buf[area][i].throughput_mbs;
	ue_ind[ue_c].x=ue_buf[area][i].x;
	ue_ind[ue_c].y=ue_buf[area][i].y;
	ue_ind[ue_c].area=area+1;
	ue_ind[ue_c].block=ue_buf[area][i].block;
	ue_ind[ue_c].block_num=ue_buf[area][i].block_num;
	ue_ind[ue_c].connect_bs_num=0;
	ue_ind[ue_c].throughput=ue_buf[area][i].throughput_mbs;
	
	for(best_c=0;best_c<cbsd_num;best_c++){
		if(best_bs<ue_buf[area][i].throughput_cbsd[best_c]){
			best_bs=ue_buf[area][i].throughput_cbsd[best_c];
			
			ue_ind[ue_c].connect_bs_num=best_c+1;
			ue_ind[ue_c].throughput=ue_buf[area][i].throughput_cbsd[best_c];
		} 
	}
	
	// add the ue connect bs number
	if(ue_ind[ue_c].connect_bs_num==0) mbs_ind.connect_num+=1;
	else cbsd_ind[(ue_ind[ue_c].connect_bs_num)-1].connect_num+=1;
	
	ue_c++;
	return ue_c;
}

//------------calculate cbsd throughput------------
float do_cbsd_pl(int bs_num, int area, int i){ // PL = 28.0 + 20log(f)|GHz| + 22log(d)|m|
	float cbsd_pl,d;
	d=sqrt(pow((float)(ue_buf[area][i].x-cbsd_ind[bs_num].x),2)+pow((float)(ue_buf[area][i].y-cbsd_ind[bs_num].y),2)); // distance unit : m
	cbsd_pl=28+10.88+22*log10(d); // in 3.5GHz
	//printf("distance : %f, pl : %f\n", d, cbsd_pl);
	return cbsd_pl;
}

float do_cbsd_P_rx(int bs_num, int area, int i, float c_p){ // Prx = Ptx + antenna gain(TX) - PL - shadow fading + antenna gain(RX)
	float P_rx,cbsd_pl;
	cbsd_pl=do_cbsd_pl(bs_num,area,i);
	if(c_p!=0){
		P_rx=(10*log10(c_p)+30)+cbsd_TX_antenna_gain-cbsd_pl-cbsd_LOS_shadow_fading+cbsd_RX_antenna_gain; // path loss unit : dB
		P_rx=pow(10.0,(P_rx-30)/10.0+9); // path loss unit : w, +9 becuz RX power is to small that can't use in float
	}else{
		P_rx=0;
	}
	//printf("cbsd power %d, Prx : %f\n", c_p, P_rx);
	return P_rx;
}

float do_cbsd_sinr(int bs_num, int area, int i, float c_p, float* c_p_r){ // SINR = Prx/(Pi+Ps)
	float SINR,P_rx,P_s,P_i=0.0;
	int P_c; // for count
	P_rx=do_cbsd_P_rx(bs_num,area,i,c_p); // unit : w
	P_s=pow(10.0,(thermal_noise-30)/10.0+9); // noise unit : w, +9 becuz noise power is to small that can't use in float

	for(P_c=0;P_c<cbsd_num;P_c++){
		if((P_c!=bs_num)&&(cbsd_ind[P_c].b_num==cbsd_ind[bs_num].b_num)) P_i+=do_cbsd_P_rx(P_c,area,i,c_p_r[P_c]);
		i++;
	}

	//printf("Prx : %f, Ps : %f\n", P_rx, P_s);
	SINR=P_rx/(P_i+P_s); // unit : w, becuz the RX and noise power and interference power have e+9, so divide it directly
	if(SINR<=0) SINR=0; // prevent C calculate fail
	//printf("Ps : %f, SNR : %f\n\n", P_s, SINR);
	return SINR;
}

float do_cbsd_throughput(int bs_num, int area, int i, float c_p, float* c_p_r, int cbsd_p_sum){
	float C,SINR;
	int bandwidth=0;

	if((c_p==0)&&(cbsd_p_sum==0)) bandwidth=0; // prevent all mbs to cbsd link are unenable
	else {
		if(cbsd_p_sum==0) cbsd_p_sum=1;
		bandwidth=(cbsd_bandwidth*cbsd_ind[bs_num].near_p_num/cbsd_p_sum)*10; // for calculate each cbsd's bandwidth
		//printf("bandwidth : %d\n",bandwidth);
		if(bandwidth==0) bandwidth=10; // bandwidth at least has 10MHz
	}
	
	SINR=do_cbsd_sinr(bs_num,area,i,c_p,&c_p_r[0]);
	C=bandwidth*(log(1+SINR)/log(2)); // unit : Mbps
	//printf("cbsd : %d, throughput: %f\n\n", bs_num, C);
	return C;
}

//------------calculate mbs throughput------------
float do_mbs_pl(int area, int i){ // LOS PL = 32.4 + 20log(f)|GHz| + 20log(d)|m| ; NLOS PL = 32.4 + 20log(f)|GHz| + 30log(d)|m|
	float mbs_pl,d;
	d=sqrt(pow((float)(ue_buf[area][i].x-mbs_ind.x),2)+pow((float)(ue_buf[area][i].y-mbs_ind.y),2)); // distance unit : m
	if(ue_buf[area][i].block!=1) mbs_pl=32.4+28.94+20*log10(d); // in 28GHz
	else mbs_pl=32.4+28.94+30*log10(d); // in 28GHz
	//printf("distance : %f, pl : %f\n", d, mbs_pl);
	return mbs_pl;
}

float do_mbs_P_rx(int area, int i, float m_p){ // Prx = Ptx + antenna gain(TX) - PL - shadow fading + antenna gain(RX)
	float P_rx,mbs_pl;
	mbs_pl=do_mbs_pl(area,i);
	if(ue_buf[area][i].block!=1) P_rx=(10*log10(m_p)+30)+mbs_TX_antenna_gain-mbs_pl-mbs_LOS_shadow_fading+mbs_RX_antenna_gain; // path loss unit : dB
	else P_rx=(10*log10(m_p)+30)+mbs_TX_antenna_gain-mbs_pl-mbs_NLOS_shadow_fading+mbs_RX_antenna_gain; // path loss unit : dB
	P_rx=pow(10.0,(P_rx-30)/10.0+9); // path loss unit : w, +9 becuz RX power is to small that can't use in float
	//printf("Prx : %f\n", P_rx);
	return P_rx;
}

float do_mbs_snr(int area, int i, float m_p){ // SNR = Prx/Ps
	float SNR,P_rx,P_s;
	P_rx=do_mbs_P_rx(area,i,m_p); // unit : w
	P_s=pow(10.0,(thermal_noise-30)/10.0+9); // noise unit : w, +9 becuz noise power is to small that can't use in float
	//printf("Prx : %f, Ps : %f\n", P_rx, P_s);
	SNR=P_rx/P_s; // unit : w, becuz both the RX and noise power have e+9, so divide it directly
	if(SNR<=0) SNR=0; // prevent C calculate fail
	//printf("Ps : %f, SNR : %f\n", P_s, SNR);
	return SNR;
}

float do_mbs_throughput(int area, int i, float m_p){ // C = B * (log_2(1+SNR))
	float C,SNR;
	SNR=do_mbs_snr(area,i,m_p);
	C=mbs_bandwidth*(log(1+SNR)/log(2)); // unit : Mbps
	//printf("throughput: %f\n\n", C);
	return C;
}

// ------------calculate throughput------------
float do_throughput(int* area_num, int area_num_max, float* m_a_p, float* c_p, float max_avg){
	float avg_t=0.0; // average throughput
	int cbsd_c=0,ue_c=0,avg_c=0; // _c for count, calculate connect cbsd throughput, max average throughput, max average times number
	int area_num_buf[4]={0}; // buffer for area number
	int i,i_c;
	int cbsd_p_sum=0; // for calculate near cbsd people number
	float cbsd_t_buf[cbsd_num]={0},cbsd_over[cbsd_num]={0}; // buffer for cbsd throughput, show the overload cbsd
									
	area_num_buf[0]=area_num[0];
	area_num_buf[1]=area_num[1];
	area_num_buf[2]=area_num[2];
	area_num_buf[3]=area_num[3];
	
	// initial ue connect bs number
	mbs_ind.connect_num=0;
	for(i=0;i<cbsd_num;i++){
		cbsd_ind[i].connect_num=0;
	}
	
	// calculate open cbsd's people
	for(i=0;i<cbsd_num;i++){
		if(c_p[i]!=0) cbsd_p_sum+=cbsd_ind[i].near_p_num;
	}

	// divide into diff area, calculate ue's throughput from connect mbs & cbsd
	for(i=0;i<area_num_max;i++){ 
		if(area_num_buf[0]!=0){
			ue_buf[0][i].throughput_mbs=do_mbs_throughput(0,i,m_a_p[0]); // connect mbs
			for(cbsd_c=0;cbsd_c<cbsd_num;cbsd_c++){
				ue_buf[0][i].throughput_cbsd[cbsd_c]=do_cbsd_throughput(cbsd_c,0,i,c_p[cbsd_c],&c_p[0],cbsd_p_sum); // connect cbsd
			}
			//ue_buf[0][i].throughput_cbsd[0]=do_cbsd_throughput(0,0,i,c1_p); // connect cbsd 1
			//ue_buf[0][i].throughput_cbsd[1]=do_cbsd_throughput(1,0,i,c2_p); // connect cbsd 2
			area_num_buf[0]-=1;
			ue_c=choose_better_bs(0,i,ue_c); // choose the better bs to connect
		}
		if(area_num_buf[1]!=0){
			ue_buf[1][i].throughput_mbs=do_mbs_throughput(1,i,m_a_p[1]);
			for(cbsd_c=0;cbsd_c<cbsd_num;cbsd_c++){
				ue_buf[1][i].throughput_cbsd[cbsd_c]=do_cbsd_throughput(cbsd_c,1,i,c_p[cbsd_c],&c_p[0],cbsd_p_sum); // connect cbsd
			}
			//ue_buf[1][i].throughput_cbsd[0]=do_cbsd_throughput(0,1,i,c1_p);
			//ue_buf[1][i].throughput_cbsd[1]=do_cbsd_throughput(1,1,i,c2_p);
			area_num_buf[1]-=1;
			ue_c=choose_better_bs(1,i,ue_c);
		}
		if(area_num_buf[2]!=0){
			ue_buf[2][i].throughput_mbs=do_mbs_throughput(2,i,m_a_p[2]);
			for(cbsd_c=0;cbsd_c<cbsd_num;cbsd_c++){
				ue_buf[2][i].throughput_cbsd[cbsd_c]=do_cbsd_throughput(cbsd_c,2,i,c_p[cbsd_c],&c_p[0],cbsd_p_sum); // connect cbsd
			}
			//ue_buf[2][i].throughput_cbsd[0]=do_cbsd_throughput(0,2,i,c1_p);
			//ue_buf[2][i].throughput_cbsd[1]=do_cbsd_throughput(1,2,i,c2_p);
			area_num_buf[2]-=1;
			ue_c=choose_better_bs(2,i,ue_c);
		}
		if(area_num_buf[3]!=0){
			ue_buf[3][i].throughput_mbs=do_mbs_throughput(3,i,m_a_p[3]);
			for(cbsd_c=0;cbsd_c<cbsd_num;cbsd_c++){
				ue_buf[3][i].throughput_cbsd[cbsd_c]=do_cbsd_throughput(cbsd_c,3,i,c_p[cbsd_c],&c_p[0],cbsd_p_sum); // connect cbsd
			}
			//ue_buf[3][i].throughput_cbsd[0]=do_cbsd_throughput(0,3,i,c1_p);
			//ue_buf[3][i].throughput_cbsd[1]=do_cbsd_throughput(1,3,i,c2_p);
			area_num_buf[3]-=1;
			ue_c=choose_better_bs(3,i,ue_c);
		}
	}

	// check whether the ue connect bs total throughput bigger than mbs to cbsd link
	for(i=0;i<cbsd_num;i++){
		cbsd_t_buf[i]=cbsd_ind[i].throughput;
	}
	// detect the overload cbsd
	for(i=0;i<ue_num;i++){
		if(ue_ind[i].connect_bs_num!=0){
			cbsd_t_buf[(ue_ind[i].connect_bs_num)-1]-=ue_ind[i].throughput;
			if(cbsd_t_buf[(ue_ind[i].connect_bs_num)-1]<0) cbsd_over[(ue_ind[i].connect_bs_num)-1]=1;
		}
	}
	// if the cbsd overload, limit the ue's throughput
	for(i=0;i<cbsd_num;i++){
		if(cbsd_over[i]==1){
			for(i_c=0;i_c<ue_num;i_c++){
				if(ue_ind[i_c].connect_bs_num==(i+1)) ue_ind[i_c].throughput=(cbsd_ind[i].throughput/cbsd_ind[i].connect_num);
			}
		}
	}
	
	// calculate average throughput
	for(i=0;i<ue_num;i++){
		avg_t=avg_t+ue_ind[i].throughput;
	}
	avg_t=avg_t/ue_num;
	//printf("a1_p : %d, a2_p : %d, a3_p : %d, a4_p : %d, c1_p : %d, c2_p : %d, average throughput : %f\n", m_a1_p, m_a2_p, m_a3_p, m_a4_p, c1_p, c2_p, avg_t);
	
	// check whether half ue bigger than average throughput
	for(i=0;i<ue_num;i++){
		if(ue_ind[i].throughput>=avg_t) avg_c+=1;
	}
	
	printf("half ue numbers %d\n\n", avg_c);
	
	//printf("half ue : %d\n", avg_c);
//	if(avg_c>=(ue_num/2)){ // if half ue bigger than average, record that time
// for test only cbsd & mbs & cbsd+mbs diff
		
		if((int)(avg_t*1000000)>(int)(max_avg*1000000)){ // record the max average power
			max_avg=avg_t;
			
			// record the bs power
			for(i=0;i<4;i++) bs_rec.a_p[i]=m_a_p[i];
			for(i=0;i<cbsd_num;i++) bs_rec.c_p[i]=c_p[i];
			bs_rec.avg_t=max_avg;
			
			// record the ue parameters
			for(i=0;i<ue_num;i++){
				ue_rec[i].x=ue_ind[i].x;
				ue_rec[i].y=ue_ind[i].y;
				ue_rec[i].area=ue_ind[i].area;
				ue_rec[i].block=ue_ind[i].block;
				ue_rec[i].block_num=ue_ind[i].block_num;
				ue_rec[i].connect_bs_num=ue_ind[i].connect_bs_num;
				ue_rec[i].throughput=ue_ind[i].throughput;
				printf("x : %3d, y : %3d, connect bs num : %d, blcok : %d, throughput : %6f\n", ue_rec[i].x, ue_rec[i].y, ue_rec[i].connect_bs_num, ue_rec[i].block, ue_rec[i].throughput);
			}
			printf("\n");
			for(i=0;i<4;i++) printf("a%d_p : %.3f, ", i+1, bs_rec.a_p[i]);
			for(i=0;i<cbsd_num;i++) printf("c%d_p : %.3f, ", i+1, bs_rec.c_p[i]);
			printf("\nmax avg : %f M bps\n\n", bs_rec.avg_t);
			printf("-----------------------------------------------\n\n");
		}
//	}
	return max_avg;
}

//------------main calculate function (allocate base station power)------------
void data_rate_calculate(int* area_num){
	float mbs_TX_power[1]={10}; // unit : w
    float cbsd_TX_power[6]={0,0.5,0.625,0.75,0.875,1}; // unit : w
	float m_a_p[4],m_c_p[cbsd_num]={0.0}; // mbs area 1~4 and mbs to cbsd power
	float c_p[cbsd_num]={0.0}; // cbsd power
	int i_m,i_c,m_c[4],c_c[cbsd_num]; // for count
	int cal_c=0;
	int area_num_max; // max ue number in 4 area
	float max_avg=0.0; // max average throughput
	float sum_p=0.0; // bs sum power
	//int r=0; // for recursive
	
	// get max area ue number in these 4 area
	if(area_num[0]<area_num[1]){
		if(area_num[1]<area_num[2]){
			if(area_num[2]<area_num[3]){
				area_num_max=area_num[3];
			}else{
				area_num_max=area_num[2];
			}
		}else{
			if(area_num[1]<area_num[3]){
				area_num_max=area_num[3];
			}else{
				area_num_max=area_num[1];
			}
		}
	}else{
		if(area_num[0]<area_num[2]){
			if(area_num[2]<area_num[3]){
				area_num_max=area_num[3];
			}else{
				area_num_max=area_num[2];
			}
		}else{
			if(area_num[0]<area_num[3]){
				area_num_max=area_num[3];
			}else{
				area_num_max=area_num[0];
			}
		}
	}
	
	// allocate bs power
	// calculate every possibility (exhaustive method)
	for(m_c[0]=0;m_c[0]<1;m_c[0]++){
		m_a_p[0]=mbs_TX_power[m_c[0]];
		for(m_c[1]=0;m_c[1]<1;m_c[1]++){
			m_a_p[1]=mbs_TX_power[m_c[1]];
			for(m_c[2]=0;m_c[2]<1;m_c[2]++){
				m_a_p[2]=mbs_TX_power[m_c[2]];
				for(m_c[3]=0;m_c[3]<1;m_c[3]++){
					m_a_p[3]=mbs_TX_power[m_c[3]];
	
					for(c_c[0]=0;c_c[0]<6;c_c[0]++){
						c_p[0]=cbsd_TX_power[c_c[0]];
						if(c_p[0]==0){
							m_c_p[0]=0; // mbs to cbsd link off
						}else{
							m_c_p[0]=mbs_cbsd_relay_p; // mbs to cbsd link on
						}
						for(c_c[1]=0;c_c[1]<6;c_c[1]++){
							c_p[1]=cbsd_TX_power[c_c[1]];
							if(c_p[1]==0){
								m_c_p[1]=0;
							}else{
								m_c_p[1]=mbs_cbsd_relay_p;
							}
							for(c_c[2]=0;c_c[2]<6;c_c[2]++){
								c_p[2]=cbsd_TX_power[c_c[2]];
								if(c_p[2]==0){
									m_c_p[2]=0;
								}else{
									m_c_p[2]=mbs_cbsd_relay_p;
								}
								for(c_c[3]=0;c_c[3]<6;c_c[3]++){
									c_p[3]=cbsd_TX_power[c_c[3]];
									if(c_p[3]==0){
										m_c_p[3]=0;
									}else{
										m_c_p[3]=mbs_cbsd_relay_p;
									}
									for(c_c[4]=0;c_c[4]<6;c_c[4]++){
										c_p[4]=cbsd_TX_power[c_c[4]];
										if(c_p[4]==0){
											m_c_p[4]=0;
										}else{
											m_c_p[4]=mbs_cbsd_relay_p;
										}
										for(c_c[5]=0;c_c[5]<6;c_c[5]++){
											c_p[5]=cbsd_TX_power[c_c[5]];
											if(c_p[5]==0){
												m_c_p[5]=0;
											}else{
												m_c_p[5]=mbs_cbsd_relay_p;
											}
											for(c_c[6]=0;c_c[6]<6;c_c[6]++){
												c_p[6]=cbsd_TX_power[c_c[6]];
												if(c_p[6]==0){
													m_c_p[6]=0;
												}else{
													m_c_p[6]=mbs_cbsd_relay_p;
												}
												for(c_c[7]=0;c_c[7]<6;c_c[7]++){
													c_p[7]=cbsd_TX_power[c_c[7]];
													if(c_p[7]==0){
														m_c_p[7]=0;
													}else{
														m_c_p[7]=mbs_cbsd_relay_p;
													}
													for(c_c[8]=0;c_c[8]<6;c_c[8]++){
														c_p[8]=cbsd_TX_power[c_c[8]];
														if(c_p[8]==0){
															m_c_p[8]=0;
														}else{
															m_c_p[8]=mbs_cbsd_relay_p;
														}
														for(c_c[9]=0;c_c[9]<6;c_c[9]++){
															c_p[9]=cbsd_TX_power[c_c[9]];
															if(c_p[9]==0){
																m_c_p[9]=0;
															}else{
																m_c_p[9]=mbs_cbsd_relay_p;
															}
																																	
															sum_p=0.0;
															printf("mbs to cbsd link check: ");
															for(i_c=0;i_c<cbsd_num;i_c++) printf("%d : %.3f, ",i_c+1,m_c_p[i_c]);
															printf("\n");
															
															for(i_m=0;i_m<4;i_m++) sum_p+=m_a_p[i_m];
															for(i_c=0;i_c<cbsd_num;i_c++) sum_p+=m_c_p[i_c];
															if(sum_p<=mbs_TX_max){ // check whether a1~a4 sum power is bigger than max power
																printf("Calculating times %d\n",cal_c++);	
																max_avg=do_throughput(&area_num[0],area_num_max,&m_a_p[0],&c_p[0],max_avg);
															}
																																								
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

//------------calculate mbs-cbsd link throughput------------
float do_mbs_cbsd_pl(int i){ // LOS PL = 32.4 + 20log(f)|GHz| + 20log(d)|km|
	float mbs_cbsd_pl,d;
	d=sqrt(pow((float)(cbsd_ind[i].x-mbs_ind.x),2)+pow((float)(cbsd_ind[i].y-mbs_ind.y),2)); // distance unit : m
	mbs_cbsd_pl=32.4+(28.94+60)+20*log10(d/1000.0); // in 28GHz
	//printf("distance : %f, pl : %f\n", d, mbs_pl);
	return mbs_cbsd_pl;
}

float do_mbs_cbsd_P_rx(int i){ // Prx = Ptx + antenna gain(TX) - PL - shadow fading + antenna gain(RX) , mbs to cbsd link Ptx = 10
	float P_rx,mbs_cbsd_pl;
	mbs_cbsd_pl=do_mbs_cbsd_pl(i);
	P_rx=(10*log10(mbs_cbsd_relay_p)+30)+mbs_TX_antenna_gain-mbs_cbsd_pl-mbs_LOS_shadow_fading+mbs_RX_antenna_gain; // path loss unit : dB
	P_rx=pow(10.0,(P_rx-30)/10.0+9); // path loss unit : w, +9 becuz RX power is to small that can't use in float
	//printf("Prx : %f\n", P_rx);
	return P_rx;
}

float do_mbs_cbsd_snr(int i){ // SNR = Prx/Ps
	float SNR,P_rx,P_s;
	P_rx=do_mbs_cbsd_P_rx(i); // unit : w
	P_s=pow(10.0,(thermal_noise-30)/10.0+9); // noise unit : w, +9 becuz noise power is to small that can't use in float
	//printf("Prx : %f, Ps : %f\n", P_rx, P_s);
	SNR=P_rx/P_s; // unit : w, becuz both the RX and noise power have e+9, so divide it directly
	if(SNR<=0) SNR=0; // prevent C calculate fail
	//printf("Ps : %f, SNR : %f\n", P_s, SNR);
	return SNR;
}

float do_mbs_cbsd_throughput(int i){ // C = B * (log_2(1+SNR))
	float C,SNR;
	SNR=do_mbs_cbsd_snr(i);
	C=mbs_cbsd_relay_bandwidth*(log(1+SNR)/log(2)); // unit : Mbps
	//printf("throughput: %f\n\n", C);
	return C;
}

//------------main deploy function------------
void deploy_node(int* area_num){
	int MAX_LINE_SIZE=40;
	char line[MAX_LINE_SIZE];
	char *result = NULL;
	int row_c=0,col_c=0,i,ue_c=0,ue_p,pepl_p; // for count
	int data[ue_num+140][4]={0}; // csv data rows , cbsd_num columns
	FILE *fpr=NULL;
	fpr=fopen("node_uniform_1200.csv","r");

    //record the data to array
    while(fgets(line, MAX_LINE_SIZE, fpr) != NULL) {
        result=strtok(line, ",");
        while(result != NULL) {
            data[row_c][col_c]=atoi(result);
            col_c++;
			result=strtok(NULL, ",");
        }
        col_c=0;
        row_c++;
    }

    // mbs coordinate
	mbs_ind.x=data[1][0];
	mbs_ind.y=data[1][1];
	printf("MBS x : %d, y : %d\n\n\n",mbs_ind.x,mbs_ind.y);

	// cbsd coordinate and channel
	for(i=0;i<cbsd_num;i++){
		cbsd_ind[i].x=data[4+i][0];
		cbsd_ind[i].y=data[4+i][1];
		cbsd_ind[i].b_num=data[4+i][2];
		cbsd_ind[i].throughput=do_mbs_cbsd_throughput(i);
		printf("CBSD %dth x : %d, y : %d, bandwidth number : %d\n\n",i,cbsd_ind[i].x,cbsd_ind[i].y,cbsd_ind[i].b_num);
	}
    printf("\n");

    // building coordinate
	for(i=0;i<(building_num/2);i++){
		building_ind[i][0].x=data[(2*i)+36][0];
		building_ind[i][0].y=data[(2*i)+36][1];
		building_ind[i][1].x=data[(2*i)+37][0];
		building_ind[i][1].y=data[(2*i)+37][1];
	}

	// area number
	for(i=0;i<4;i++){
		area_num[i]=data[30+30*2+8][i];
		printf("area number %d : %d\n\n",i+1,area_num[i]);
	}
	printf("\n");

	// ue parameter
	ue_p=30+30*2+12;
	for(i=0;i<4;i++){

		while(ue_c<area_num[i]){
			ue_buf[i][ue_c].x=data[ue_p+ue_c][0];
			ue_buf[i][ue_c].y=data[ue_p+ue_c][1];
			ue_buf[i][ue_c].block=data[ue_p+ue_c][2];
			ue_buf[i][ue_c].block_num=data[ue_p+ue_c][3];
			printf("ue x : %d, y : %d, block : %d, block cbsd number : %d\n\n", ue_buf[i][ue_c].x, ue_buf[i][ue_c].y, ue_buf[i][ue_c].block, ue_buf[i][ue_c].block_num);
			ue_c++;
		}
		ue_c=0;
		ue_p=ue_p+area_num[i]+1;
	}
	printf("\n");

	// near cbsd people number
	pepl_p=30+30*2+ue_num+17;
	for(i=0;i<cbsd_num;i++){
		cbsd_ind[i].near_p_num=data[pepl_p+i][0]; // near cbsd people number
	}

	//check cbsd people number
	printf("\ncbsd people number\n");
	for(i=0;i<cbsd_num;i++){
		printf("%d: %d,",i+1,cbsd_ind[i].near_p_num);
	}
	printf("\n");

    fclose (fpr);
}

int main(){
	int area_num[4]={0};
	
	deploy_node(&area_num[0]); // deploy every node
	
	data_rate_calculate(&area_num[0]); // calculate ue's data rate
	
	write_excel(); // write the result into excel
	
	printf("hello world\n");
	
	return 0;
}
