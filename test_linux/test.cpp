#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <algorithm>
#include "daggen_commons.h"
#include "mcts.h"



using namespace std;

static DAG generateDAG(void);
static void generateTasks(DAG dag);
static void generateDependencies(DAG dag);
static void generateTransfers(DAG dag);
static void freeDAG(DAG dag);
static void init_comm_costs(DAG dag);
static void init_user_satisfaction();

double avg_mcts_makespan = 0;
double avg_mcts_satisfaction = 0;
//double avg_mcts_speedup = 0;
//double avg_mcts_runtime = 0;

double avg_ccr = 0;
double avg_comm_cost = 0;
double avg_comp_cost = 0;

int times = 500;
double start, finish;
FILE* fp_result;
int main(int argc, char** argv)
{

    srand((unsigned int)getpid() + (unsigned int)time(NULL));
    /* parse command line options */
    if (parseOptions(argc, argv) == -1)
    {
        printUsage();
        exit(1);
    }
    
    fp_result = fopen("fp_result.txt", "w");
    fprintf(fp_result, "Iteration_numbers:%d\n", global.n_simulations);
    fprintf(fp_result, "User_numbers:%d, Subtask_numbers:%d, Dependency_degree:%d\n", global.nb_users, global.nb_tasks_per_user, global.nb_dependency);
    fprintf(fp_result, "Server_numbers:%d\n", global.n_processors);
    fprintf(fp_result, "Data_size:%d - %d, Comp_to_Comm_Cycle:%d\n", MIN_cost, MAX_cost,  Comp_to_Comm_Cycle);
    
    for (int kk = 0; kk < times; kk++)
    {
        printf("times = %d\n", kk + 1);
        DAG dag;
        Queue queue;
        Tree tree;
  
        init_processor();

        init_bandwidth();

        assign_processor_bandwidth();      

        dag = generateDAG();

        init_user_satisfaction();

        queue = init_queue();
        tree = init_tree(dag, queue);
        start = clock();
        UCTsearch(dag, tree, queue);
        finish = clock();



        printf( "CCR:%lf, comm_cost_per_task:%lf, comp_cost_per_task:%lf\n", global.avg_comm_latency_per_task / global.avg_comp_latency_per_task, global.avg_comm_latency_per_task, global.avg_comp_latency_per_task);
        printf( "mcts_makespan = %lf, mcts_satisfaction = %lf", result.mcts_makespan, result.mcts_satisfaction);
        printf("\n");

        avg_mcts_makespan += result.mcts_makespan;
        avg_mcts_satisfaction += result.mcts_satisfaction;
       


        avg_ccr += global.avg_comm_latency_per_task / global.avg_comp_latency_per_task;
        avg_comm_cost += global.avg_comm_latency_per_task;
        avg_comp_cost += global.avg_comp_latency_per_task;

        freeDAG(dag);
        //free(tree);
        //free(queue);
        free(global.user_satisfaction);
        free(mcts_g.processor_performance);
        for (int i = 0; i < global.n_processors; i++)
           free( mcts_g.bandwidth[i]);
        free(mcts_g.bandwidth);
        for (int i = 0; i < global.n * global.n; i++)
            free(mcts_g.comm_costs[i]);
        free(mcts_g.comm_costs);
        free(mcts_g.nodetrace);
        for (int i = 0; i < global.n + 1; i++)
            free(mcts_g.scheduled[i]);
        free(mcts_g.scheduled);
        for (int i = 0; i < global.n_processors; i++)
            free(mcts_g.scheduled_sort_on_per_server[i]);
        free(mcts_g.scheduled_sort_on_per_server);
     //   fprintf(fp_result, "\n");
    }
      fprintf(fp_result, "avg_ccr = %lf,avg_comm_cost = %lf, avg_comp_cost = %lf, avg_cost= %lf\n", avg_ccr / times, avg_comm_cost / times, avg_comp_cost / times, avg_comm_cost / times + avg_comp_cost / times);
      fprintf(fp_result, "avg_mcts_makespan = %lf, avg_mcts_satisfaction = %lf", avg_mcts_makespan / times, avg_mcts_satisfaction / times);
      fprintf(fp_result, "\n\n");
      printf("\n\n");

      printf("avg_ccr = %lf,avg_comm_cost = %lf, avg_comp_cost = %lf, avg_cost= %lf\n", avg_ccr / times, avg_comm_cost / times, avg_comp_cost / times, avg_comm_cost / times + avg_comp_cost / times);
      printf( "avg_mcts_makespan = %lf, avg_mcts_satisfaction = %lf\n", avg_mcts_makespan / times, avg_mcts_satisfaction / times);
     
      fclose(fp_result);
      exit(0);
}

void init_user_satisfaction()
{
    global.user_satisfaction = (double*)calloc(global.nb_users, sizeof(double));
    //double avg_satisfaction = (global.avg_comm_latency_per_task * (global.nb_tasks_per_user - 1) + global.avg_comp_latency_per_task * global.nb_tasks_per_user) / global.n_processors;
    double avg_satisfaction = 1;
    for (int i = 0; i <= global.nb_users - 1; i++)
    {
        global.user_satisfaction[i] = 7.0 * rand() / RAND_MAX +1.0;
    }
    sort(global.user_satisfaction, global.user_satisfaction + global.nb_users);
    printf("user_satisfaction:");
    for (int i = 0; i < global.nb_users; i++)
        printf("%lf ", global.user_satisfaction[i]);
    printf("\n");
}


static DAG generateDAG(void)
{

    DAG dag;

    dag = (DAG)calloc(1, sizeof(struct _DAG));
    
    /* Generating all the tasks */
    generateTasks(dag);

    /* Generating the Dependencies */
    generateDependencies(dag);

    /* Generating the transfer costs */
    generateTransfers(dag);

    //computeDownRank(dag);

    //computeUpRank(dag);

    //computeSumRank(dag);

    //computeCP(dag);

    //generateCP(dag);


    return dag;
}

static void generateTasks(DAG dag)
{
    
    int nb_levels= global.nb_tasks_per_user - global.nb_dependency + 1;
    int* nb_tasks = NULL;
    nb_tasks = (int*)realloc(nb_tasks, (nb_levels + 1) * sizeof(int));   
   
    global.n = 0;

    for (int i = 0; i < nb_levels; i++)
    {
        if(i== (nb_levels - 2))
            nb_tasks[i] = global.nb_users * global.nb_dependency;
        else
            nb_tasks[i] = global.nb_users;
        global.n += nb_tasks[i];
    }

    dag->nb_levels = nb_levels;
    dag->levels = (Task**)calloc(dag->nb_levels, sizeof(Task*));
    dag->nb_tasks_per_level = nb_tasks;
    int tmp_priority = MAX_priority-1;
    for (int i = 0; i < dag->nb_levels; i++)
    {
        dag->levels[i] = (Task*)calloc(dag->nb_tasks_per_level[i], sizeof(Task));       
        for (int j = 0; j < dag->nb_tasks_per_level[i]; j++)
        {
            dag->levels[i][j] = (Task)calloc(1, sizeof(struct _Task));
            dag->levels[i][j]->cost = getRandomNumberBetween(MIN_cost, MAX_cost);
            dag->levels[i][j]->data_size = dag->levels[i][j]->cost / Comp_to_Comm_Cycle;
            if (i == nb_levels - 2)
                dag->levels[i][j]->priority = tmp_priority - j / global.nb_dependency;
            else
            dag->levels[i][j]->priority = tmp_priority-j;
        }      
    }

}

static void generateDependencies(DAG dag)
{
    int nb_parents;
    Task parent;
    for (int i = 1; i < dag->nb_levels; i++)
    {
        for (int j = 0; j < dag->nb_tasks_per_level[i]; j++)
        {
            if (i == dag->nb_levels - 2)
            {
                nb_parents = 1;
                dag->levels[i][j]->nb_parents = nb_parents;
                dag->levels[i][j]->nb_parents_r = nb_parents;
                dag->levels[i][j]->parents = (Task*)calloc(nb_parents, sizeof(Task));
                parent = dag->levels[i - 1][(int)(j / global.nb_dependency)];
                dag->levels[i][j]->parents[0] = parent;   
                parent->children = (Task*)realloc(parent->children, (parent->nb_children + 1) * sizeof(Task));
                parent->children[(parent->nb_children)] = dag->levels[i][j];
          
                (parent->nb_children)++;
            }
            else if (i == dag->nb_levels - 1)
            {
                nb_parents = global.nb_dependency;
                dag->levels[i][j]->nb_parents = nb_parents;
                dag->levels[i][j]->nb_parents_r = nb_parents;
                dag->levels[i][j]->parents = (Task*)calloc(nb_parents, sizeof(Task));
                for (int k = 0; k < nb_parents; k++)
                {
                    parent = dag->levels[i - 1][j * global.nb_dependency + k];
                    dag->levels[i][j]->parents[k] = parent;
                    parent->children = (Task*)realloc(parent->children, (parent->nb_children + 1) * sizeof(Task));
                    parent->children[(parent->nb_children)] = dag->levels[i][j];
                    (parent->nb_children)++;
                }
                
            }
            else
            {
                nb_parents = 1;
                dag->levels[i][j]->nb_parents = nb_parents;
                dag->levels[i][j]->nb_parents_r = nb_parents;
                dag->levels[i][j]->parents = (Task*)calloc(nb_parents, sizeof(Task));
                parent = dag->levels[i - 1][j];
                dag->levels[i][j]->parents[0] = parent;
                parent->children = (Task*)realloc(parent->children, (parent->nb_children + 1) * sizeof(Task));
                parent->children[(parent->nb_children)] = dag->levels[i][j];
                (parent->nb_children)++;
            }
        }
    }
    for (int i = 0; i < dag->nb_levels; i++)
    {
        for (int j = 0; j < dag->nb_tasks_per_level[i]; j++)
        {
            dag->levels[i][j]->comm_costs = (double*)calloc(
                dag->levels[i][j]->nb_children, sizeof(double));
            dag->levels[i][j]->transfer_tags = (int*)calloc(
                dag->levels[i][j]->nb_children, sizeof(int));
            dag->levels[i][j]->comp_costs = (double*)calloc(
                global.n_processors, sizeof(double));
            //dag->levels[i][j]->MAB_rewards = (double*)calloc(
            //    global.n_processors, sizeof(double));
            //dag->levels[i][j]->MAB_numbers = (int*)calloc(
            //    global.n_processors, sizeof(int));
        }
    }
    int k = 0;
    global.nb_parents = (int*)calloc(global.n* global.n, sizeof(int));
    for (int i = 1; i < dag->nb_levels; i++)
        for (int j = 0; j < dag->nb_tasks_per_level[i]; j++)
        {
            global.nb_parents[k] = dag->levels[i][j]->nb_parents;
            k++;
        }
    int node_count = 1;
    for (int i = 0; i < dag->nb_levels; i++)
    {
        for (int j = 0; j < dag->nb_tasks_per_level[i]; j++)
        {
            dag->levels[i][j]->tag = node_count;
            node_count++;
        }
    }

}

static void generateTransfers(DAG dag)// Enforces the CCR
{
    int i, j, k;
    double temp;
    double sum_comp_cost = 0;
    int times = 0;
    for (i = 0; i < dag->nb_levels; i++)
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
        {
            temp = 0;
            for (k = 0; k < global.n_processors; k++)
            {
                dag->levels[i][j]->comp_costs[k] =
                    dag->levels[i][j]->cost / mcts_g.processor_performance[k];
                temp += dag->levels[i][j]->comp_costs[k];
                
            }
            dag->levels[i][j]->ave_comp_cost = temp / global.n_processors;
            sum_comp_cost += dag->levels[i][j]->cost / mcts_g.avg_processor_performance;
            times++;
        }
    global.avg_comp_latency_per_task = sum_comp_cost;
    for (i = 0; i < dag->nb_levels - 1; i++)
    {
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
        {
            for (k = 0; k < dag->levels[i][j]->nb_children; k++)
            {
                dag->levels[i][j]->comm_costs[k] = dag->levels[i][j]->data_size;
            }
        }
    }
     init_comm_costs(dag); 
    return;
}



void init_comm_costs(DAG dag) {
    int i, j, k, p_id, c_id;
    double sum_comm_cost = 0;
    int times = 0;
    for (i = 0; i < dag->nb_levels - 1; i++)
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
            for (k = 0; k < dag->levels[i][j]->nb_children; k++)
            {
                Task task = dag->levels[i][j];
                p_id = task->tag;
                c_id = task->children[k]->tag;
                mcts_g.comm_costs[p_id][c_id] = task->comm_costs[k];
                sum_comm_cost += 8.0 * task->comm_costs[k]/((global.n_MBS_processors* MBS_trans_rate + (global.n- global.n_MBS_processors)* mcts_g.avg_bandwidth*log2(1 + global.P_BS * global.g_bs / double(global.Noise)))/global.n);
                times++;
            }
        }
    global.avg_comm_latency_per_task = sum_comm_cost;


}


void freeDAG(DAG dag)
{
    int i, j;
    for (i = 0; i < dag->nb_levels; i++)
    {
        for (j = 0; j < dag->nb_tasks_per_level[i]; j++)
        {
            free(dag->levels[i][j]->transfer_tags);
            free(dag->levels[i][j]->children);
            free(dag->levels[i][j]->comm_costs);
            free(dag->levels[i][j]);
        }
        free(dag->levels[i]);
    }
    free(dag->levels);
    free(dag->nb_tasks_per_level);
    free(dag);
}