#pragma once
#ifndef MCTS_H_
#define MCTS_H_

//#define MINP 900 //
//#define MAXP 1000 //
#define Comp_to_Comm_Cycle 200
#define MIN_cost 900//
#define MAX_cost 1000

#define User_CPU_frequency 1000
#define MIN_SBS_CPU_frequency 4000
#define MAX_SBS_CPU_frequency 5000
#define MIN_MBS_CPU_frequency 7000
#define MAX_MBS_CPU_frequency 8000

#define User_bandwidth 5
#define MIN_SBS_bandwidth 20
#define MAX_SBS_bandwidth 25
#define MBS_trans_rate 100.0 //wired transmission rate
#define MAXtask 500
#define MAX_priority 10000
//#define fp_process (mcts_g.output_process)
typedef struct _Node* Node;
typedef struct _Tree* Tree;
typedef struct _Queue* Queue;
typedef struct _Element* Element;
typedef struct _scheduled* Scheduled;
//typedef struct _TempScheduled *TempScheduled;

typedef struct {
    double* processor_performance;
    double** bandwidth;
    double avg_processor_performance;
    double avg_bandwidth;
    double** comm_costs;
    double makespan;
    int node_id;
    int nb_node;
    Node* nodetrace;
    Scheduled* scheduled;
    //double* availP;
    int count;
    int* ans;
    int* visit_mark;
    Task* schedule_Task;
    int** edge_mark;
    int** scheduled_sort_on_per_server;
    //int count; 
    //FILE *output_process;
//    int **comp_costs;
//    int **comm_costs;
} MCTS_g;

extern MCTS_g mcts_g;

struct _scheduled {
    Task task;
    int processor;
    int priority;
    double EST;//start time
    double EFT;//finish time
};


struct _Element {
    Task task;
    Element next;
};

struct _Queue {
    int n;
    Element head;
};


struct _Node {
    int node_id;
    int task_id;
    Task task;
    int processor;
    double Q;//reward
    int visits;// visit numbers
    int nb_children;
    Node* children;
    int marker;
};


struct _Tree {
    int depth;
    Node root;
};
void print_mcts_info(Tree tree, Queue queue);
void expand(Node node, Queue queue);
void assign_processor_bandwidth();
double get_comm_time(int parent_id, int parent_proc, int task_id, int processor);
//void Assign_processor_bandwidth();
//void init_comm_costs(DAG dag);
Task init_root_task(DAG dag);
void UCTsearch(DAG dag, Tree tree, Queue queue);
Tree init_tree(DAG dag, Queue queue);
Queue init_queue();
void init_processor();
//void BB_init_processor();
void init_bandwidth();
//void BB_init_processor();
void refresh_time_on_succedtask(int task_id, double add_time);
void refresh_time_on_processor(int task_id, int processor);
void update_scheduled_task(Task task, int processor, FILE* fp);
double get_makespan();
int get_best_processor(Task task);
double get_avg_satisfaction(DAG dag);
Task delete_n_queue(Queue queue, int n);
Task delete_task_queue(Queue queue, int task_id);
void insert_queue(Queue queue, Task task);
void update_queue(Queue queue, Task* task_ready, int nb_tasks);
void init_children(Node node, Queue queue);
Node best_child(Node node);
void update_ready_task(DAG dag, Task task, Queue queue);
double get_start_time(Task task, int processor);
double find_start_time(Task task, int processor);
void update_scheduled_sort_on_processor(int task_id);
void update_scheduled_task(Task task, int processor, FILE* fp);
void backup();
Node tree_policy(DAG dag, Tree tree, Queue queue);
void default_random_policy(DAG dag, Node leafnode, Queue queue);
void output_tree(Tree tree);
void output_subtree(Node node, FILE* fp);
void print_queue(Queue queue);
Tree init_tree(DAG dag, Queue queue);
double get_speedup(DAG dag, FILE* fp);
#endif /*MCTS_H_*/
