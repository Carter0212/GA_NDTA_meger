#define pi 3.1415926

#define ue_num 1200
#define cbsd_num 30
#define building_num 30
#define implement_Time 49

// increase mbs_TX_max value for extend topology
#define mbs_TX_max 75 // unit w
#define mbs_ue_max 10 // unit w
#define mbs_TX_antenna_gain 20 // unit dbi
#define mbs_RX_antenna_gain 9 // unit dbi
#define mbs_LOS_shadow_fading 4 // unit db
#define mbs_NLOS_shadow_fading 7.8 // unit db
#define mbs_bandwidth 100
#define mbs_cbsd_relay_p 7.5 // unit w
#define mbs_cbsd_relay_bandwidth 1200

#define cbsd_TX_antenna_gain 5 // unit dbi
#define cbsd_RX_antenna_gain 5 // unit dbi
#define cbsd_LOS_shadow_fading 3 // unit db
#define cbsd_bandwidth 15

#define thermal_noise -83.02

#define population_size 300
#define crossover_prob 85 // means crossover probability = 85%
#define mutation_prob 1 // means mutation probability = 1%
#define iteration_num 100 // max generation number
#define elitism_num 5 // the elitism number

// for strategy II
#define mbs_crossover_prob 90 // means mbs part crossover probability
#define cbsd_crossover_prob 0 // means cbsd part crossover probability

#define init_temp 500 // SA Initial temperature
#define cooling_fac 0.8 // SA Cooling Factor

// for water filling
#define wf_threshold_l 7.24 // water filling lower threshold
#define wf_threshold_u 100 // water filling upper threshold

//------------parameters define------------

struct GA_max_avg_throughput_generatopm_record{
	float number; // resord the max avg. throughput number (generation)
	int number_buf; // buffer for resord the max avg. throughput number (generation) 
	float throughput; // record the max avg. throughput
};

struct chromosome_buffer{
	int area[4][5]; // area 1~4 power choose
	int mbs_cbsd[cbsd_num]; // mbs to cbsd connect, 0=disconnect, 1=connect 
	int cbsd[cbsd_num][5]; // cbsd power choose
	int number; // resord chromosome number
};

struct chromosome_parameters{
	int area[4][5]; // area 1~4 power choose
	int mbs_cbsd[cbsd_num]; // mbs to cbsd connect, 0=disconnect, 1=connect 
	int cbsd[cbsd_num][5]; // cbsd power choose
	int number; // resord chromosome number
};

struct first_ue_record{
	int x; // ue's x axis
	int y; // ue's y axis
	int area; // which mbs area
	int block; // judge whether NLOS with mbs, 1=NLOS , 0=LOS
	int block_num; // be blocked cbsd number
	int connect_bs_num; // connect bs number, 0=mbs , 1=cbsd1 , 2=cbsd2
	float throughput;
};

struct max_ue_record{
	int x; // ue's x axis
	int y; // ue's y axis
	int area; // which mbs area
	int block; // judge whether NLOS with mbs, 1=NLOS , 0=LOS
	int block_num; // be blocked cbsd number
	int connect_bs_num; // connect bs number, 0=mbs , 1=cbsd1 , 2=cbsd2
	float throughput;
};

struct max_bs_record{
	float a_p[4]; // area number
	float c_p[cbsd_num];
	float avg_t; // average throughput
	int times[5]; // calculate times in each thread
};

struct ue_buffer{
	int x; // ue's x axis
	int y; // ue's y axis
	int block; // judge whether NLOS with mbs, 1=NLOS , 0=LOS
	int block_num; // be blocked cbsd number
	float throughput_mbs; // connect mbs throughtput
	float throughput_cbsd[cbsd_num]; // connect cbsd throughtput
};

struct ue_parameters{
	int x; // ue's x axis
	int y; // ue's y axis
	int area; // which mbs area
	int block; // judge whether NLOS with mbs, 1=NLOS , 0=LOS
	int block_num; // be blocked cbsd number
	int connect_bs_num; // connect bs number, 0=mbs , 1=cbsd1 , 2=cbsd2 ...
	float throughput; // ue's throughput
};

struct cbsd_parameters{
	int x; // cbsd's x axis
	int y; // cbsd's y axis
	int b_num; // cbsd's bandwidth number
	int connect_num; // connect each cbsd number
	int near_p_num; // near cbsd people number
	float throughput; // mbs to cbsd link throughput
	int enable_channel[cbsd_bandwidth];//cbsd use which one bandwith total have 15 channel
};

struct building_parameters{
	int x; // building's x axis
	int y; // building's y axis
};

struct mbs_parameters{
	int x; // mbs's x axis
	int y; // mbs's y axis
	int connect_num; // connect mbs number
	int cbsd_same_ch_num; // record the number of cbsd which is same channel
};

struct bandwidth_parameters{
	int Channel_Usage_Situation[cbsd_bandwidth]; // us to update & cheack Remaining channels
	
};
