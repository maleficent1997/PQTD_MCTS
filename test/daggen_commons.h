/******************************************************************************
 * Copyright (c) 2007-2013. F. Suter, S. Hunold.
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL 2.1) which comes with this package.
 *****************************************************************************/

 //#include <cstdio>
#ifndef DAGGEN_COMMONS_H_
#define DAGGEN_COMMONS_H_

//#ifndef MAX
//#  define MAX(x,y) ((x) > (y) ? (x) : (y))
//#endif
//#ifndef MIN
//#  define MIN(x,y) ((x) < (y) ? (x) : (y))
//#endif



#define OUTPUT (global.output_file)
#define GIGA 1024*1024*1024

/*********************************/
/** Global variables            **/
/*********************************/

typedef struct {
	int n;          /* number of tasks in the graph */
	int n_processors; /* number of processors in the system*/
	int n_MBS_processors;/* n_processors中 MBS处理器的数量*/
	int n_simulations;
	int cp;//UCT公式中衡量Exploitation 和 Exploration的参数Cpuct
	double fat;     /* fatness parameter                    */
	double regular; /* regularity                         */
	double ccr;     /* Communication to computation ratio */
	double density;
	double mindata, maxdata;
	double minalpha, maxalpha; /* Amdahl's law parameter */
	int jump;
	int dot_output;
	int* nb_parents;
	FILE* output_file;
	int nb_users;//用户数量
	int nb_tasks_per_user;
	int nb_dependency;
	double distance;
	double User_BS;
	double P_BS;
	double g_bs;
	double Noise;
	double* user_satisfaction;
	double avg_comp_latency_per_task;
	double avg_comm_latency_per_task;
} Global;
extern Global global;

typedef struct _Task* Task;
typedef struct _DAG* DAG;

struct _Task {
	int tag;
	double cost;
	double* comp_costs;
	double data_size;
	int priority;
	double alpha;
	int nb_children;
	int nb_parents;
	int nb_parents_r;
	Task* children;
	Task* parents;
	double* comm_costs;
	int* transfer_tags;
	double urank;
	double drank;
	double sum_rank;
	double cp;
	double ave_comp_cost;
	int flat;
};

struct _DAG {
	int nb_levels;//level numbers
	int* nb_tasks_per_level;//subtask numbers per level
	Task** levels;//[nb_levels][nb_tasks_per_level]
};

typedef struct {

	double mcts_makespan;
	double mcts_satisfaction;
	double mcts_runtime;
	double mcts_speedup;

}Result;

extern Result result;

int parseOptions(int argc, char* const* argv);


void printUsage(void);

void outputDAG(DAG dag);
void outputDOT(DAG dag);
void print_dag_info(DAG dag);

int getIntRandomNumberAround(int x, double perc);

double getRandomNumberBetween(double x, double y);
//double getRandomNumberBetween(double x, double y) {
//	double r;
//
//	r = x + (y - x) * rand() / (RAND_MAX + 1.0);
//	return r;
//}

void generate_randomID(int min, int max, int n, int a[]);

void IDtoIJ(int id, int* ij, DAG dag);

int get_nb_tasks(DAG dag, int n_level);

// void saveDAG(DAG dag);

// DAG readDAG();

#endif /*DAGGEN_COMMONS_H_*/

