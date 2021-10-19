#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "node.h"
#define deploy_flag 0 // if true, uniform user distribution , else, non-uniform user distribution
#define block_flag 0 // if true, real block judge , else, define a block probability
#define block_prob 25

struct ue_buffer ue_buf[4][ue_num];
struct ue_parameters ue_ind[ue_num];
struct cbsd_parameters cbsd_ind[building_num];
struct building_parameters building_ind[building_num][2];
struct mbs_parameters mbs_ind;

//------------calculate near cbsd people number------------
void peoeple_cal(){
	int i,j;
	float ue_bs_dis=0.0;
	for(i=0;i<cbsd_num;i++){
		cbsd_ind[i].near_p_num=0; // near cbsd people number
	}
	for(i=0;i<ue_num;i++){
		for(j=0;j<cbsd_num;j++){
			ue_bs_dis=pow((pow(ue_ind[i].x-cbsd_ind[j].x,2)+pow(ue_ind[i].y-cbsd_ind[j].y,2)),0.5);
			if(ue_bs_dis<=500) cbsd_ind[j].near_p_num++;
		}
	}

	printf("near cbsd people number : ");
	for(i=0;i<cbsd_num;i++){
		printf("%d = %d, ",i+1,cbsd_ind[i].near_p_num);
	}
	printf("\n\n");
}


//------------judge ue's area------------
int judge_area(int i){
	if(ue_ind[i].x>=0) {
		if(ue_ind[i].y>=0){
			return 1; // area 1
		}else{
			return 4; // area 4
		}
	} else {
		if(ue_ind[i].y>=0){
			return 2; // area 2
		}else{
			return 3; // area 3
		}
	}
}

//------------judge block------------
int judge_block(int i){ // equation : y=ax+b
	float b_a[building_num],b_b[building_num]; // building a,b
	float u_a,u_b; // ue with MBS a,b
	float x[building_num],y[building_num]; // cross between building and UE
	int ab_c,xy_c;

	// calculate the building equation
	for(ab_c=0;ab_c<building_num;ab_c++){
		b_a[ab_c]=(float)(building_ind[ab_c][0].y-building_ind[ab_c][1].y)/(float)(building_ind[ab_c][0].x-building_ind[ab_c][1].x);
		b_b[ab_c]=(float)(building_ind[ab_c][0].y*building_ind[ab_c][1].x-building_ind[ab_c][1].y*building_ind[ab_c][0].x)/(float)(building_ind[ab_c][1].x-building_ind[ab_c][0].x);
		//printf("cbsd %d, a : %f, b : %f\n", ab_c, b_a[ab_c], b_b[ab_c]);
	}

	// ue to mbs equation
	u_a=(float)(ue_ind[i].y-mbs_ind.y)/(float)(ue_ind[i].x-mbs_ind.x);
	u_b=(float)(ue_ind[i].y*mbs_ind.x-mbs_ind.y*ue_ind[i].x)/(float)(mbs_ind.x-ue_ind[i].x);
	//printf("ue's ab , a : %f, b : %f\n", u_a, u_b);

	// find the cross point in these two equations
	for(xy_c=0;xy_c<building_num;xy_c++){
		x[xy_c]=(b_b[xy_c]-u_b)/(u_a-b_a[xy_c]);
		y[xy_c]=(u_b*b_a[xy_c]-b_b[xy_c]*u_a)/(b_a[xy_c]-u_a);
		//printf("%d line corss , x : %f, y : %f\n", xy_c, x[xy_c], y[xy_c]);

		// check whether the cross point is in the interval (between ue and mbs && building)
		if((((y[xy_c]-mbs_ind.y)*(y[xy_c]-ue_ind[i].y))<=0.0) && ((x[xy_c]-building_ind[xy_c][1].x)*(x[xy_c]-building_ind[xy_c][0].x))<=0.0) {
			ue_ind[i].block_num=xy_c+1; // record the block building number
			return 1; // ue be blocked
		}
		//printf("?? 1: %f, 2 : %f\n", (y[xy_c]-mbs_ind.y)*(y[xy_c]-ue_ind[i].y), (x[xy_c]-building_ind[xy_c][1].x)*(x[xy_c]-building_ind[xy_c][0].x));
	}
	ue_ind[i].block_num=0;
	return 0; // ue not be blocked
}

//------------deploy ue------------
void ue_deploy(int* area_num){
	int i,block_n=0;
	int a1_c=0,a2_c=0,a3_c=0,a4_c=0;
	int deploy_f=0,deploy_c;
	float ue_bs_dis=0.0;
	srand(time(NULL));
	for(i=0;i<ue_num;i++){ // deploy ue
		if(deploy_flag==1) { // uniform distribute
			ue_ind[i].x=(rand() % 4000)-2000;
			ue_ind[i].y=(rand() % 4000)-2000;
			if((ue_ind[i].x==0)&&(ue_ind[i].y==0)){ // prevent same as mbs
				ue_ind[i].x=(rand() % 4000)-2000;
				ue_ind[i].y=(rand() % 4000)-2000;
			}
		} else { // non-uniform distribute
			while(deploy_f==0) {
				ue_ind[i].x=(rand() % 4000)-2000;
				ue_ind[i].y=(rand() % 4000)-2000;
				ue_bs_dis=pow((pow(ue_ind[i].x-mbs_ind.x,2)+pow(ue_ind[i].y-mbs_ind.y,2)),0.5);
				if(ue_bs_dis<200.0){
					if((rand()%100)<((-2*pow(ue_bs_dis,2)/pow(200,2)+1)*100)) deploy_f=1; // P=((-2 * d(x,y)^2) / D(x,y)^2) + 1 , d means the distance between ue and bs , D means the max bs coverage
				}
				if(deploy_f==0){
					for(deploy_c=0;deploy_c<cbsd_num;deploy_c++) {
						ue_bs_dis=pow((pow(ue_ind[i].x-cbsd_ind[deploy_c].x,2)+pow(ue_ind[i].y-cbsd_ind[deploy_c].y,2)),0.5);
						//printf("distance : %f\n", ue_bs_dis);
						if(ue_bs_dis<1500.0){
							if((rand()%100)<((-2*pow(ue_bs_dis,2)/pow(1500,2)+1)*100)) deploy_f=1; // P=((-2 * d(x,y)^2) / D(x,y)^2) + 1 , d means the distance between ue and bs , D means the max bs coverage
						}
					}
				}
			}
			deploy_f=0;
		}
		
		ue_ind[i].area=judge_area(i); // judge ue's area
		
		if(block_flag==1) {
			ue_ind[i].block=judge_block(i); // judge whether ue be blocked by building
		} else {
			if((rand()%100)<block_prob) ue_ind[i].block=1;
		}
		
		if(ue_ind[i].block==1) block_n+=1;
		printf("ue x : %d, y : %d, area : %d, block : %d, block cbsd number : %d\n\n\n", ue_ind[i].x, ue_ind[i].y, ue_ind[i].area, ue_ind[i].block, ue_ind[i].block_num);

		// use a buffer structure to divide ue into each area
		if(ue_ind[i].area==1) {
			ue_buf[0][a1_c].x=ue_ind[i].x;
			ue_buf[0][a1_c].y=ue_ind[i].y;
			ue_buf[0][a1_c].block=ue_ind[i].block;
			ue_buf[0][a1_c].block_num=ue_ind[i].block_num;
			a1_c+=1;
		} else if(ue_ind[i].area==2) {
			ue_buf[1][a2_c].x=ue_ind[i].x;
			ue_buf[1][a2_c].y=ue_ind[i].y;
			ue_buf[1][a2_c].block=ue_ind[i].block;
			ue_buf[1][a2_c].block_num=ue_ind[i].block_num;
			a2_c+=1;
		} else if(ue_ind[i].area==3) {
			ue_buf[2][a3_c].x=ue_ind[i].x;
			ue_buf[2][a3_c].y=ue_ind[i].y;
			ue_buf[2][a3_c].block=ue_ind[i].block;
			ue_buf[2][a3_c].block_num=ue_ind[i].block_num;
			a3_c+=1;
		} else if(ue_ind[i].area==4) {
			ue_buf[3][a4_c].x=ue_ind[i].x;
			ue_buf[3][a4_c].y=ue_ind[i].y;
			ue_buf[3][a4_c].block=ue_ind[i].block;
			ue_buf[3][a4_c].block_num=ue_ind[i].block_num;
			a4_c+=1;
		}
	}

	// record each area ue number
	area_num[0]=a1_c;
	area_num[1]=a2_c;
	area_num[2]=a3_c;
	area_num[3]=a4_c;

	printf("---------------------------\n\n\n");
	printf("area number 1: %d, 2: %d, 3: %d, 4: %d", area_num[0], area_num[1], area_num[2], area_num[3]);
	printf("\nblock ue number : %d", block_n);
	printf("\n\n\n---------------------------\n\n\n");
}

//------------deploy building------------
void building_deploy(){
	building_ind[0][0].x=cbsd_ind[0].x+85;
	building_ind[0][0].y=cbsd_ind[0].y+85;
	building_ind[0][1].x=cbsd_ind[0].x-85;
	building_ind[0][1].y=cbsd_ind[0].y-85;
	building_ind[1][0].x=cbsd_ind[1].x+90;
	building_ind[1][0].y=cbsd_ind[1].y-90;
	building_ind[1][1].x=cbsd_ind[1].x-90;
	building_ind[1][1].y=cbsd_ind[1].y+90;
	building_ind[2][0].x=cbsd_ind[2].x+80;
	building_ind[2][0].y=cbsd_ind[2].y;
	building_ind[2][1].x=cbsd_ind[2].x-80;
	building_ind[2][1].y=cbsd_ind[2].y;
	building_ind[3][0].x=cbsd_ind[3].x-80;
	building_ind[3][0].y=cbsd_ind[3].y+80;
	building_ind[3][1].x=cbsd_ind[3].x+80;
	building_ind[3][1].y=cbsd_ind[3].y-80;
	building_ind[4][0].x=cbsd_ind[4].x+82;
	building_ind[4][0].y=cbsd_ind[4].y+70;
	building_ind[4][1].x=cbsd_ind[4].x-82;
	building_ind[4][1].y=cbsd_ind[4].y-70;
	building_ind[5][0].x=cbsd_ind[5].x-85;
	building_ind[5][0].y=cbsd_ind[5].y+75;
	building_ind[5][1].x=cbsd_ind[5].x+85;
	building_ind[5][1].y=cbsd_ind[5].y-75;
	building_ind[6][0].x=cbsd_ind[6].x+72;
	building_ind[6][0].y=cbsd_ind[6].y+76;
	building_ind[6][1].x=cbsd_ind[6].x-72;
	building_ind[6][1].y=cbsd_ind[6].y-76;
	building_ind[7][0].x=cbsd_ind[7].x-95;
	building_ind[7][0].y=cbsd_ind[7].y+65;
	building_ind[7][1].x=cbsd_ind[7].x+95;
	building_ind[7][1].y=cbsd_ind[7].y-65;
	building_ind[8][0].x=cbsd_ind[8].x-55;
	building_ind[8][0].y=cbsd_ind[8].y-85;
	building_ind[8][1].x=cbsd_ind[8].x+55;
	building_ind[8][1].y=cbsd_ind[8].y+85;
	building_ind[9][0].x=cbsd_ind[9].x-63;
	building_ind[9][0].y=cbsd_ind[9].y+61;
	building_ind[9][1].x=cbsd_ind[9].x+63;
	building_ind[9][1].y=cbsd_ind[9].y-61;
	building_ind[10][0].x=cbsd_ind[10].x+98;
	building_ind[10][0].y=cbsd_ind[10].y+94;
	building_ind[10][1].x=cbsd_ind[10].x-98;
	building_ind[10][1].y=cbsd_ind[10].y-94;
	building_ind[11][0].x=cbsd_ind[11].x-75;
	building_ind[11][0].y=cbsd_ind[11].y-90;
	building_ind[11][1].x=cbsd_ind[11].x+75;
	building_ind[11][1].y=cbsd_ind[11].y+90;
	building_ind[12][0].x=cbsd_ind[12].x-98;
	building_ind[12][0].y=cbsd_ind[12].y-47;
	building_ind[12][1].x=cbsd_ind[12].x+98;
	building_ind[12][1].y=cbsd_ind[12].y+47;
	building_ind[13][0].x=cbsd_ind[13].x+90;
	building_ind[13][0].y=cbsd_ind[13].y-49;
	building_ind[13][1].x=cbsd_ind[13].x-90;
	building_ind[13][1].y=cbsd_ind[13].y+49;
	building_ind[14][0].x=cbsd_ind[14].x+75;
	building_ind[14][0].y=cbsd_ind[14].y+75;
	building_ind[14][1].x=cbsd_ind[14].x-75;
	building_ind[14][1].y=cbsd_ind[14].y-75;
	building_ind[15][0].x=cbsd_ind[15].x-63;
	building_ind[15][0].y=cbsd_ind[15].y+98;
	building_ind[15][1].x=cbsd_ind[15].x+63;
	building_ind[15][1].y=cbsd_ind[15].y-98;
	building_ind[16][0].x=cbsd_ind[16].x+60;
	building_ind[16][0].y=cbsd_ind[16].y+60;
	building_ind[16][1].x=cbsd_ind[16].x-60;
	building_ind[16][1].y=cbsd_ind[16].y-60;
	building_ind[17][0].x=cbsd_ind[17].x-90;
	building_ind[17][0].y=cbsd_ind[17].y+90;
	building_ind[17][1].x=cbsd_ind[17].x+90;
	building_ind[17][1].y=cbsd_ind[17].y-90;
	building_ind[18][0].x=cbsd_ind[18].x+75;
	building_ind[18][0].y=cbsd_ind[18].y+75;
	building_ind[18][1].x=cbsd_ind[18].x-75;
	building_ind[18][1].y=cbsd_ind[18].y-75;
	building_ind[19][0].x=cbsd_ind[19].x-77;
	building_ind[19][0].y=cbsd_ind[19].y+77;
	building_ind[19][1].x=cbsd_ind[19].x+77;
	building_ind[19][1].y=cbsd_ind[19].y-77;
	building_ind[20][0].x=cbsd_ind[20].x+50;
	building_ind[20][0].y=cbsd_ind[20].y+90;
	building_ind[20][1].x=cbsd_ind[20].x-60;
	building_ind[20][1].y=cbsd_ind[20].y-90;
	building_ind[21][0].x=cbsd_ind[21].x-45;
	building_ind[21][0].y=cbsd_ind[21].y+85;
	building_ind[21][1].x=cbsd_ind[21].x+45;
	building_ind[21][1].y=cbsd_ind[21].y-85;
	building_ind[22][0].x=cbsd_ind[22].x+43;
	building_ind[22][0].y=cbsd_ind[22].y+83;
	building_ind[22][1].x=cbsd_ind[22].x-43;
	building_ind[22][1].y=cbsd_ind[22].y-83;
	building_ind[23][0].x=cbsd_ind[23].x-88;
	building_ind[23][0].y=cbsd_ind[23].y+45;
	building_ind[23][1].x=cbsd_ind[23].x+88;
	building_ind[23][1].y=cbsd_ind[23].y-45;
	building_ind[24][0].x=cbsd_ind[24].x+100;
	building_ind[24][0].y=cbsd_ind[24].y+46;
	building_ind[24][1].x=cbsd_ind[24].x-100;
	building_ind[24][1].y=cbsd_ind[24].y-46;
	building_ind[25][0].x=cbsd_ind[25].x-100;
	building_ind[25][0].y=cbsd_ind[25].y+100;
	building_ind[25][1].x=cbsd_ind[25].x+100;
	building_ind[25][1].y=cbsd_ind[25].y-100;
	building_ind[26][0].x=cbsd_ind[26].x+100;
	building_ind[26][0].y=cbsd_ind[26].y+100;
	building_ind[26][1].x=cbsd_ind[26].x-100;
	building_ind[26][1].y=cbsd_ind[26].y-100;
	building_ind[27][0].x=cbsd_ind[27].x-100;
	building_ind[27][0].y=cbsd_ind[27].y+100;
	building_ind[27][1].x=cbsd_ind[27].x+100;
	building_ind[27][1].y=cbsd_ind[27].y-100;
	building_ind[28][0].x=cbsd_ind[28].x+100;
	building_ind[28][0].y=cbsd_ind[28].y+100;
	building_ind[28][1].x=cbsd_ind[28].x-100;
	building_ind[28][1].y=cbsd_ind[28].y-100;
	building_ind[29][0].x=cbsd_ind[29].x-100;
	building_ind[29][0].y=cbsd_ind[29].y+100;
	building_ind[29][1].x=cbsd_ind[29].x+100;
	building_ind[29][1].y=cbsd_ind[29].y-100;
}

//------------deploy cbsd------------
void cbsd_deploy(){
	cbsd_ind[0].x=-970;
	cbsd_ind[0].y=900;
	cbsd_ind[0].b_num=2; // bandwidth number
	cbsd_ind[1].x=740;
	cbsd_ind[1].y=900;
	cbsd_ind[1].b_num=1;
	cbsd_ind[2].x=980;
	cbsd_ind[2].y=-860;
	cbsd_ind[2].b_num=1;
	cbsd_ind[3].x=-750;
	cbsd_ind[3].y=-850;
	cbsd_ind[3].b_num=5;
	cbsd_ind[4].x=-750;
	cbsd_ind[4].y=1000;
	cbsd_ind[4].b_num=4;
	cbsd_ind[5].x=1050;
	cbsd_ind[5].y=600;
	cbsd_ind[5].b_num=2;
	cbsd_ind[6].x=760;
	cbsd_ind[6].y=-1135;
	cbsd_ind[6].b_num=3;
	cbsd_ind[7].x=-870;
	cbsd_ind[7].y=-1100;
	cbsd_ind[7].b_num=1;
	cbsd_ind[8].x=-1080;
	cbsd_ind[8].y=450;
	cbsd_ind[8].b_num=1;
	cbsd_ind[9].x=970;
	cbsd_ind[9].y=530;
	cbsd_ind[9].b_num=3;
	cbsd_ind[10].x=540;
	cbsd_ind[10].y=-1000;
	cbsd_ind[10].b_num=2;
	cbsd_ind[11].x=-1200;
	cbsd_ind[11].y=-300;
	cbsd_ind[11].b_num=3;
	cbsd_ind[12].x=-200;
	cbsd_ind[12].y=1170;
	cbsd_ind[12].b_num=7;
	cbsd_ind[13].x=820;
	cbsd_ind[13].y=1255;
	cbsd_ind[13].b_num=4;
	cbsd_ind[14].x=930;
	cbsd_ind[14].y=-1080;
	cbsd_ind[14].b_num=6;
	cbsd_ind[15].x=-1115;
	cbsd_ind[15].y=-505;
	cbsd_ind[15].b_num=2;
	cbsd_ind[16].x=-180;
	cbsd_ind[16].y=1230;
	cbsd_ind[16].b_num=3;
	cbsd_ind[17].x=500;
	cbsd_ind[17].y=1150;
	cbsd_ind[17].b_num=8;
	cbsd_ind[18].x=1200;
	cbsd_ind[18].y=-450;
	cbsd_ind[18].b_num=10;
	cbsd_ind[19].x=-500;
	cbsd_ind[19].y=-1100;
	cbsd_ind[19].b_num=7;
	cbsd_ind[20].x=-650;
	cbsd_ind[20].y=980;
	cbsd_ind[20].b_num=9;
	cbsd_ind[21].x=1180;
	cbsd_ind[21].y=450;
	cbsd_ind[21].b_num=9;
	cbsd_ind[22].x=700;
	cbsd_ind[22].y=-1160;
	cbsd_ind[22].b_num=7;
	cbsd_ind[23].x=-1050;
	cbsd_ind[23].y=-900;
	cbsd_ind[23].b_num=8;
	cbsd_ind[24].x=-1200;
	cbsd_ind[24].y=820;
	cbsd_ind[24].b_num=10;
	cbsd_ind[25].x = 1200;
    cbsd_ind[25].y = 1000;
	cbsd_ind[25].b_num=11;
    cbsd_ind[26].x = 1340;
    cbsd_ind[26].y = -1120;
	cbsd_ind[26].b_num=12;
    cbsd_ind[27].x = -1440;
    cbsd_ind[27].y = -1280;
	cbsd_ind[27].b_num=10;
    cbsd_ind[28].x = -1160;
    cbsd_ind[28].y = 1350;
	cbsd_ind[28].b_num=13;
    cbsd_ind[29].x = 1420;
    cbsd_ind[29].y = 960;
	cbsd_ind[29].b_num=14;
}

//------------deploy mbs------------
void mbs_deploy(){
	mbs_ind.x=0;
	mbs_ind.y=0;
}

//------------main deploy function------------
void deploy_node(int* area_num,FILE *fp){
	int i,a_c=0; // for count
	mbs_deploy(); // deploy mbs
	cbsd_deploy(); // deploy cbsd (in the middle of the building)
	building_deploy(); // deploy the building
	ue_deploy(&area_num[0]); // deploy ue and fill the struct
	peoeple_cal(); // calculate near cbsd people number

	fprintf(fp,"mbs coordinate\n");
	fprintf(fp,"%d,%d\n\n",mbs_ind.x,mbs_ind.y);
	fprintf(fp,"cbsd coordinate and channel\n");
	for(i=0;i<cbsd_num;i++){
		fprintf(fp,"%d,%d,%d\n",cbsd_ind[i].x,cbsd_ind[i].y,cbsd_ind[i].b_num);
	}
	fprintf(fp,"\n");
	fprintf(fp,"building coordinate\n");
	for(i=0;i<building_num;i++){
		fprintf(fp,"%d,%d\n",building_ind[i][0].x,building_ind[i][0].y);
		fprintf(fp,"%d,%d\n",building_ind[i][1].x,building_ind[i][1].y);
	}
	fprintf(fp,"\n");
	fprintf(fp,"area numbers\n");
	fprintf(fp,"%d,%d,%d,%d\n\n",area_num[0],area_num[1],area_num[2],area_num[3]);
	fprintf(fp,"x , y , block , block number\n");
	for(i=0;i<4;i++){
		fprintf(fp,"area %d\n",i+1);
		while(a_c<area_num[i]){
			fprintf(fp,"%d,%d,%d,%d\n",ue_buf[i][a_c].x,ue_buf[i][a_c].y,ue_buf[i][a_c].block,ue_buf[i][a_c].block_num);
			a_c++;
		}
		a_c=0;
	}

	fprintf(fp,"\nnear cbsd peopel number\n");
	for(i=0;i<cbsd_num;i++){
		fprintf(fp,"%d\n",cbsd_ind[i].near_p_num);
	}

}


int main(){
	int area_num[4]={0};
	int time_count;
	char *filename = NULL;
	FILE *fp=NULL;
	char data; //data all mix
	char data_name_unifom[] = "node_uniform_1200_";
	char data_name_nonuniform[] = "node_nonuniform_1200_";
	char data_type[] = ".csv";
	char data_num;
	for(time_count=0; time_count<implement_Time ; time_count++){
		if(deploy_flag==1) {
			sprintf(data_num, "%d", time_count);
        	strcpy(data,data_name_unifom);
        	strcat(data,buf);
        	strcat(data,data_type);
        	filename = data;
        	printf("now start make date %s\n", filename);
			fp=fopen("node_uniform_1200.csv","w");
		} 
		else {
			sprintf(data_num, "%d", time_count);
        	strcpy(data,data_name_nonuniform);
        	strcat(data,buf);
        	strcat(data,data_type);
        	filename = data;
        	printf("now start make date %s\n", filename);
			fp=fopen("node_nonuniform_1200.csv","w");
		}

	deploy_node(&area_num[0],fp,time_count); // deploy every node
	fclose(fp);

	}

	


	printf("hello world\n");

	return 0;
}
