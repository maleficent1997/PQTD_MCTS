/******************************************************************************
 * Copyright (c) 2007-2013. F. Suter, S. Hunold.
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL 2.1) which comes with this package.
 *****************************************************************************/
#pragma once;
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
 //#include "getopt.h"
#include "daggen_commons.h"

 /*********************************/
 /** Command line option parsing **/
 /*********************************/

 /*
  * Parse the options and set default values
  *
  * returns -1 for usage
  * returns -2 for no usage
  * returns 0 when Ok
  */

Global global;
Result result;
int parseOptions(int argc, char* const* argv) {
    int ret_val = 0;
    int c;

    int oflag = 0;
    int nflag = 0;
    int pflag = 0;
    int sflag = 0;
    int cflag = 0;
    int fat_flag = 0;
    int density_flag = 0;
    int ccr_flag = 0;
    int mindata_flag = 0;
    int maxdata_flag = 0;
    int minalpha_flag = 0;
    int maxalpha_flag = 0;
    int regular_flag = 0;
    int jump_flag = 0;
    int dot_flag = 0;
    FILE* fp;
    fp = fopen("digraph.dot", "w");
    global.output_file = fp;
    global.jump = 1;
    global.mindata = 2048;
    global.maxdata = 11264;
    global.minalpha = 0.0;
    global.maxalpha = 0.2;
    global.ccr = 0.2;
    global.density = 0.5;
    global.fat = 0.5;
    global.regular = 0.9;

    global.n_processors = 5;
    global.n_MBS_processors = 2;
    global.n_simulations = 1000;
    global.cp = -10;
    global.dot_output = 0;

    global.nb_users = 5;
    global.nb_tasks_per_user = 10;
    global.nb_dependency = 5;
    global.n = global.nb_tasks_per_user * global.nb_users;
    global.distance = 50.0;
    global.User_BS = 0.2; 
    global.P_BS = 20.0; 
    global.g_bs = 127 + 30 * log(global.distance); 
    global.Noise = 1e-13; 
    return ret_val;
}


/*
 * printUsage()
 */
void printUsage(void)
{
    fprintf(stderr, "daggen [options] [-o <output file>]\n"
        "\t -n <number of tasks>\n"
        "\t --mindata <minimum data size>\n"
        "\t --maxdata <maximum data size>\n"
        "\t --minalpha <minimum Amdahl's law parameter value>\n"
        "\t --maxalpha <maximum Amdahl's law parameter value>\n");
    fprintf(stderr,
        "\t --fat <dag shape>\n"
        "\t	 fat=1.0: fat (maximum parallelism)\n"
        "\t	 fat=0.0: thin (minimum parallelism)\n"
        "\t --density <density>\n"
        "\t    density=0.0: minimum number of dependencies\n"
        "\t    density=1.0: full graph\n"
        "\t --regular <regularity for num tasks per level>\n"
        "\t    regular= 1.0: perfectly regular\n"
        "\t    regular= 0.0: irregular\n"
        "\t --ccr <communication(MBytes) to computation(sec) ratio>\n"
        "\t --jump <number of levels spanned by communications>\n"
        "\t    jump=1: perfectly synchronized levels\n"
        "\t --dot: output generated DAG in the DOT format\n"
        "\n");
    return;
}


/****************/
/** DAG output **/
/****************/

void outputDAG(DAG dag) {
    int i, j, k;
    /* starting at 1 for the root node */
    int node_count = 1;

    /* count and tag the nodes */
    for (i = 0; i < dag->nb_levels; i++) {
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            dag->levels[i][j]->tag = node_count++;
        }
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            for (k = 0; k < dag->levels[i][j]->nb_children; k++) {
                dag->levels[i][j]->transfer_tags[k] = node_count++;
            }
        }
    }
    /* accounting for the END node */
    fprintf(OUTPUT, "NODE_COUNT %d\n", node_count + 1);

    /* Create the root node */
    fprintf(OUTPUT, "NODE 0 ");
    for (i = 0; i < dag->nb_tasks_per_level[0] - 1; i++) {
        fprintf(OUTPUT, "%d,", dag->levels[0][i]->tag);
    }
    if (dag->nb_tasks_per_level[0])
        fprintf(OUTPUT, "%d ROOT 0.0 0.0\n", dag->levels[0][i]->tag);
    else
        fprintf(OUTPUT, "%d ROOT 0.0 0.0\n", node_count);

    /* Creating the regular nodes until next to last level */
    for (i = 0; i < dag->nb_levels - 1; i++) {
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            /* do the COMPUTATION */
            fprintf(OUTPUT, "NODE %d ", dag->levels[i][j]->tag);
            for (k = 0; k < dag->levels[i][j]->nb_children - 1; k++) {
                fprintf(OUTPUT, "%d,", dag->levels[i][j]->transfer_tags[k]);
            }
            if (dag->levels[i][j]->nb_children) {
                fprintf(OUTPUT, "%d COMPUTATION %.0f %.2f\n",
                    dag->levels[i][j]->transfer_tags[k],
                    dag->levels[i][j]->cost,
                    dag->levels[i][j]->alpha);
            }
            else {
                fprintf(OUTPUT, "%d COMPUTATION %.0f %.2f\n",
                    node_count,
                    dag->levels[i][j]->cost,
                    dag->levels[i][j]->alpha);
            }
            /* do the TRANSFER */
            for (k = 0; k < dag->levels[i][j]->nb_children; k++) {
                fprintf(OUTPUT, "NODE %d ", dag->levels[i][j]->transfer_tags[k]);
                fprintf(OUTPUT, "%d TRANSFER %.0f 0.0\n",
                    dag->levels[i][j]->children[k]->tag,
                    dag->levels[i][j]->comm_costs[k]);
            }
        }
    }

    /* Do the last level */
    for (j = 0; j < dag->nb_tasks_per_level[dag->nb_levels - 1]; j++) {
        fprintf(OUTPUT, "NODE %d %d COMPUTATION %.0f %.2f\n",
            dag->levels[dag->nb_levels - 1][j]->tag,
            node_count,
            dag->levels[dag->nb_levels - 1][j]->cost,
            dag->levels[dag->nb_levels - 1][j]->alpha);
    }

    /* Do the end node */
    fprintf(OUTPUT, "NODE %d - END 0.0 0.0\n", node_count);
}

void outputDOT(DAG dag) {
    int i, j, k;
    /* starting at 1 for the root node */
    int node_count = 1;

    /* count and tag the nodes */
    for (i = 0; i < dag->nb_levels; i++) {
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            dag->levels[i][j]->tag = node_count++;
        }
        //    for (j=0; j<dag->nb_tasks_per_level[i]; j++) {
        //      for (k=0; k<dag->levels[i][j]->nb_children; k++) {
        //        dag->levels[i][j]->transfer_tags[k] = node_count++;
        //      }
        //    }
    }
    /* accounting for the END node */
    fprintf(OUTPUT, "digraph G {\n");

    /* Create the root node */
    fprintf(OUTPUT, "  0 [label=root,size=\"0\",alpha=\"0\"]\n");
    for (i = 0; i < dag->nb_tasks_per_level[0]; i++) {
        fprintf(OUTPUT, "  0 -> %d [label = \"0\"]\n", dag->levels[0][i]->tag);
    }

    /* Creating the regular nodes until next to last level */
    for (i = 0; i < dag->nb_levels; i++) {
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            /* do the COMPUTATION */
            fprintf(OUTPUT, "  %d [label= \"%d C=%.3lf\",size=\"%.0f\", alpha=\"%.2f\"]\n",
                dag->levels[i][j]->tag,
                dag->levels[i][j]->tag - 1,
                dag->levels[i][j]->cost,
                dag->levels[i][j]->cost,
                dag->levels[i][j]->alpha);

            /* do the TRANSFER */
            for (k = 0; k < dag->levels[i][j]->nb_children; k++) {
                fprintf(OUTPUT, "  %d -> %d [label =\"%.f\"]\n",
                    dag->levels[i][j]->tag,
                    dag->levels[i][j]->children[k]->tag,
                    dag->levels[i][j]->comm_costs[k]);
            }
        }
    }

    //  /* Do the last level */
    //  for (j=0; j<dag->nb_tasks_per_level[dag->nb_levels-1]; j++) {
    //    fprintf(OUTPUT,"  %d -> %d [label = \"0\"]\n",
    //        dag->levels[dag->nb_levels-1][j]->tag,
    //        node_count);
    //  }

    for (i = 0; i < dag->nb_levels; i++)
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            if (dag->levels[i][j]->nb_children == 0) {
                //   fprintf(OUTPUT, "0 child£º%d\n",dag->levels[i][j]->tag);
                fprintf(OUTPUT, "  %d -> %d [label = \"0\"]\n",
                    dag->levels[i][j]->tag, node_count);
            }
        }
    fprintf(OUTPUT, "}\n");
    fclose(OUTPUT);
}

void print_dag_info(DAG dag) {
    int i, j, k;

    printf("task costs:\n");
    for (i = 0; i < dag->nb_levels; i++) {
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            printf("%.2f ", dag->levels[i][j]->cost);
        }
        printf("\n");
    }
    printf("task comp_costs:\n");
    for (i = 0; i < dag->nb_levels; i++)
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            for (k = 0; k < global.n_processors; k++) {
                printf("%.2f ", dag->levels[i][j]->comp_costs[k]);
            }
            printf("\n");
        }
    printf("task datasize:\n");
    for (i = 0; i < dag->nb_levels; i++) {
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            printf("%lf ", dag->levels[i][j]->data_size);
        }
        printf("\n");
    }

    printf("task comm_costs:\n");
    for (i = 0; i < dag->nb_levels; i++) {
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            for (k = 0; k < dag->levels[i][j]->nb_children; k++) {
                printf("%.2f ", dag->levels[i][j]->comm_costs[k]);
            }
            printf("\n");
        }
    }
}

/***********************/
/** Random generators **/
/***********************/

/*
 * getIntRandomNumberAround()
 *
 * returns a STRICTLY POSITIVIE int around x, by perc percents.
 */
int getIntRandomNumberAround(int x, double perc) {
    double r;
    int new_int;
    r = -perc + (2 * perc * rand() / (RAND_MAX + 1.0));
    int kk = (int)((double)x * (1.0 + r / 100.00));
    new_int = (1 > kk) ? 1 : kk;
    return new_int;
}

/*
 * getRandomNumberBetween()
 *
 */
double getRandomNumberBetween(double x, double y) {
    double r;

    r = x + (y - x) * rand() / (RAND_MAX + 1.0);
    return r;
}

// generate n number of non-repetitive random IDs between min and max
void generate_randomID(int min, int max, int n, int a[]) {

    srand(time(NULL));

    int RandNum, i, j, flag = 0, t = 0;

    while (1) {
        flag = 0;
        if (t == n)
            break;
        RandNum = (rand() % (max - min)) + min;
        for (i = 0; i < t; i++) {
            if (a[i] == RandNum)
                flag = 1;
        }
        if (flag != 1) {
            a[t] = RandNum;
            // printf("%d ",a[t]);
            t++;
        }
    }
    //printf("\n");
}
// given node ID, get level i, index j
void IDtoIJ(int id, int* ij, DAG dag) {
    int i, j;
    for (i = 0; i < dag->nb_levels; i++)
    {
        if (id < dag->nb_tasks_per_level[i])
        {
            break;
        }
        else
        {
            id = id - dag->nb_tasks_per_level[i];
        }
    }
    ij[0] = i;
    ij[1] = id;
    //  printf("ij[]:%d,%d\n",ij[0],ij[1]);
}

// void saveDAG(DAG dag){
//   FILE *fp,*fp1;
//   int g_size;
//   fp = fopen("dag_struct.ds","w");
//   fp1 = fopen("dag_size.txt","w");
//   g_size = sizeof(dag);
//   fwrite(&g_size,sizeof(int),1,fp1);
//   fwrite(dag, sizeof(dag),1,fp);
//   fclose(fp);
// }

// DAG readDAG(){
//   FILE *fp,*fp1;
//   DAG *dag;
//   int g_size;
//   fp1 = fopen("dag_size.txt","r");
//   fp = fopen("dag_struct.ds","r");
//   fread(&g_size,sizeof(int),1,fp1);
//   printf("g_size:%d\n",g_size);
//   dag = (DAG)calloc(1,sizeof(g_size));
//   fread(dag, g_size, 1, fp);
//   return dag;
// }

int get_nb_tasks(DAG dag, int n_level) {
    int i;
    int n_tasks = 0;
    for (i = 0; i < n_level; i++) {
        n_tasks += dag->nb_tasks_per_level[i];
    }
    return n_tasks;
}

void set_computation_cost(DAG dag, int n_proc, int W_dag, double beta) {
    int i, j;
    int w_i;
    srand(time(NULL));
    for (i = 0; i < dag->nb_levels; i++)
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {

        }
}

