#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <cstdio>
#include "daggen_commons.h"
#include "mcts.h"
#define MAXtime 50
#define CP -10
#define _CRT_SECURE_NO_WARNINGS
MCTS_g mcts_g;


void init_processor()
{
    int i;
    FILE* fp_pc;
    errno_t err_kk = fopen_s(&fp_pc, "processor.txt", "w");

    if (err_kk != 0) {
        fprintf(stderr, "Error: Unable to open file 'processor.txt'.\n");
        return;
    }

    // Print the User CPU frequency
    fprintf(fp_pc, "%lf ", (double)User_CPU_frequency);

    // Print SBS CPU frequencies
    for (i = 1; i < global.n_processors - global.n_MBS_processors; i++)
    {
        fprintf(fp_pc, "%lf ", getRandomNumberBetween((double)MIN_SBS_CPU_frequency, (double)MAX_SBS_CPU_frequency));
    }

    // Print MBS CPU frequencies
    for (i = global.n_processors - global.n_MBS_processors; i < global.n_processors; i++)
    {
        fprintf(fp_pc, "%lf ", getRandomNumberBetween((double)MIN_MBS_CPU_frequency, (double)MAX_MBS_CPU_frequency));
    }

    fprintf(fp_pc, "\n");
    fclose(fp_pc);
}


void init_bandwidth() {
    int i, j;
    FILE* fp_bw;
    errno_t err = fopen_s(&fp_bw, "bandwidth.txt", "w");
    //fp_bw = fopen("bandwidth.txt", "w");
    for (i = 0; i < global.n_processors; i++) {
        for (j = 0; j < global.n_processors; j++) {
            if (i != j)
            {
                if (j == 0)
                    fprintf(fp_bw, "%lf ", (double)User_bandwidth);
                else
                    fprintf(fp_bw, "%lf ", getRandomNumberBetween((double)MIN_SBS_bandwidth, (double)MAX_SBS_bandwidth));
            }
            else {
                fprintf(fp_bw, "0 ");
            }
        }
        fprintf(fp_bw, "\n");
    }
    fprintf(fp_bw, "\n");
    fclose(fp_bw);
}


void assign_processor_bandwidth()
{
    int i, j;
    FILE* fp_pc, * fp_bw;
    double sum_processor_performance = 0;
    double sum_bandwidth = 0;
    errno_t err = fopen_s(&fp_pc, "processor.txt", "r");
    //fp_pc = fopen("processor.txt", "r");
    errno_t err_ps = fopen_s(&fp_bw, "bandwidth.txt", "r");
    //fp_bw = fopen("bandwidth.txt", "r");
    if (!fp_pc || !fp_bw)
    {
        mcts_g.processor_performance = (double*)calloc(global.n_processors, sizeof(double));
        mcts_g.bandwidth = (double**)calloc(global.n_processors, sizeof(double*));
        for (i = 0; i < global.n_processors; i++)
            mcts_g.bandwidth[i] = (double*)calloc(global.n_processors, sizeof(double));
        // mcts_g.comm_costs = (double**)calloc(global.n, sizeof(double*));
        mcts_g.comm_costs = (double**)calloc(global.n * global.n, sizeof(double*));
        for (i = 0; i < global.n * global.n; i++)
            // mcts_g.comm_costs[i] = (double*)calloc(global.n, sizeof(double));
            mcts_g.comm_costs[i] = (double*)calloc(global.n * global.n, sizeof(double));
        mcts_g.makespan = 0; mcts_g.nb_node = 0; mcts_g.node_id = 0;
        mcts_g.nodetrace = (Node*)calloc(global.n + 1, sizeof(Node));
        // mcts_g.scheduled = (Scheduled*)calloc(global.n + 1, sizeof(Scheduled));
        mcts_g.scheduled = (Scheduled*)calloc(global.n + 1, sizeof(Scheduled));
        // mcts_g.availP = (double*)calloc(global.n_processors, sizeof(double));
        for (i = 0; i < global.n + 1; i++)
            mcts_g.scheduled[i] = (Scheduled)calloc(1, sizeof(struct _scheduled));
        mcts_g.scheduled_sort_on_per_server = (int**)calloc(global.n_processors, sizeof(int*));
        for (i = 0; i < global.n_processors; i++)
            mcts_g.scheduled_sort_on_per_server[i] = (int*)calloc(global.n + 1, sizeof(int));
        mcts_g.processor_performance[0] = (double)User_CPU_frequency;

        for (i = 1; i < global.n_processors - global.n_MBS_processors; i++)
            mcts_g.processor_performance[i] = getRandomNumberBetween((double)MIN_SBS_CPU_frequency, (double)MAX_SBS_CPU_frequency);
        for (i = global.n_processors - global.n_MBS_processors; i < global.n_processors; i++)
            mcts_g.processor_performance[i] = getRandomNumberBetween((double)MIN_MBS_CPU_frequency, (double)MAX_MBS_CPU_frequency);
        for (i = 0; i < global.n_processors; i++)
            sum_processor_performance += mcts_g.processor_performance[i];
        mcts_g.avg_processor_performance = sum_processor_performance / global.n_processors;
        for (i = 0; i < global.n_processors; i++) {
            for (j = 0; j < global.n_processors; j++)
            {
                if (i != j)
                {
                    if (j == 0)
                        mcts_g.bandwidth[i][j] = (double)User_bandwidth;
                    else
                        mcts_g.bandwidth[i][j] = getRandomNumberBetween((double)MIN_SBS_bandwidth, (double)MAX_SBS_bandwidth);
                }
                else {
                    mcts_g.bandwidth[i][j] = 0;
                }
                sum_bandwidth += mcts_g.bandwidth[i][j];
            }

        }
        mcts_g.avg_bandwidth = sum_bandwidth / ((global.n_processors - global.n_MBS_processors) * (global.n_processors - global.n_MBS_processors));
    }
    else
    {
        //Allocate memory space
        mcts_g.processor_performance = (double*)calloc(global.n_processors, sizeof(double));
        mcts_g.bandwidth = (double**)calloc(global.n_processors, sizeof(double*));
        for (i = 0; i < global.n_processors; i++)
            mcts_g.bandwidth[i] = (double*)calloc(global.n_processors, sizeof(double));
        // mcts_g.comm_costs = (double**)calloc(global.n, sizeof(double*));
        mcts_g.comm_costs = (double**)calloc(global.n * global.n, sizeof(double*));
        for (i = 0; i < global.n * global.n; i++)
            // mcts_g.comm_costs[i] = (double*)calloc(global.n, sizeof(double));
            mcts_g.comm_costs[i] = (double*)calloc(global.n * global.n, sizeof(double));
        mcts_g.makespan = 0; mcts_g.nb_node = 0; mcts_g.node_id = 0;
        mcts_g.nodetrace = (Node*)calloc(global.n + 1, sizeof(Node));
        // mcts_g.scheduled = (Scheduled*)calloc(global.n + 1, sizeof(Scheduled));
        mcts_g.scheduled = (Scheduled*)calloc(global.n + 1, sizeof(Scheduled));
        // mcts_g.availP = (double*)calloc(global.n_processors, sizeof(double));
        for (i = 0; i < global.n + 1; i++)
            mcts_g.scheduled[i] = (Scheduled)calloc(1, sizeof(struct _scheduled));
        mcts_g.scheduled_sort_on_per_server = (int**)calloc(global.n_processors, sizeof(int*));
        for (i = 0; i < global.n_processors; i++)
            mcts_g.scheduled_sort_on_per_server[i] = (int*)calloc(global.n + 1, sizeof(int));
        //initialize mcts_g.scheduled_sort_on_per_server
        //for (i = 0; i < global.n_processors; i++)
        //    for (j = 0; j < global.n_processors; j++)
        //        mcts_g.scheduled_sort_on_per_server[i][j] = INT_MAX;

        for (i = 0; i < global.n_processors; i++)
        {
            fscanf_s(fp_pc, "%lf", &mcts_g.processor_performance[i]);
            //fscanf(fp_pc, "%lf ", &mcts_g.processor_performance[i]);
            sum_processor_performance += mcts_g.processor_performance[i];
            //printf("%lf ", mcts_g.processor_performance[i]);
        }
        mcts_g.avg_processor_performance = sum_processor_performance / global.n_processors;

        // printf("\n");

        for (i = global.n_MBS_processors; i < global.n_processors; i++) {
            for (j = global.n_MBS_processors; j < global.n_processors; j++) {
                fscanf_s(fp_bw, "%lf ", &mcts_g.bandwidth[i][j]);
                //fscanf(fp_bw, "%lf ", &mcts_g.bandwidth[i][j]);
                sum_bandwidth += mcts_g.bandwidth[i][j];
            }
            mcts_g.avg_bandwidth = sum_bandwidth / ((global.n_processors - global.n_MBS_processors) * (global.n_processors - global.n_MBS_processors));
            fscanf_s(fp_bw, "\n");
            //fscanf(fp_bw, "\n");
        }
    }
}

double get_comm_time(int parent_id, int parent_proc, int task_id, int processor)
{
    double comm_time = 0;
    if (mcts_g.bandwidth[parent_proc][processor] != 0) {
        if (parent_proc >= 0 && parent_proc < global.n_processors - global.n_MBS_processors && processor >= 0 && processor < global.n_processors - global.n_MBS_processors)
        {
            if (parent_proc == 0)
                comm_time = 8.0 * mcts_g.comm_costs[parent_id][task_id] / (mcts_g.bandwidth[parent_proc][processor] * log2(1 + global.User_BS * global.g_bs / double(global.Noise)));
            else
                comm_time = 8.0 * mcts_g.comm_costs[parent_id][task_id] / (mcts_g.bandwidth[parent_proc][processor] * log2(1 + global.P_BS * global.g_bs / double(global.Noise)));
        }
        else if (parent_proc == 0 && processor >= global.n_processors - global.n_MBS_processors && processor < global.n_processors)
            comm_time = 8.0 * mcts_g.comm_costs[parent_id][task_id] / (mcts_g.bandwidth[parent_proc][processor] * log2(1 + global.User_BS * global.g_bs / double(global.Noise))) + 8.0 * mcts_g.comm_costs[parent_id][task_id] / MBS_trans_rate;
        else if (processor == 0 && parent_proc >= global.n_processors - global.n_MBS_processors && parent_proc < global.n_processors)
            comm_time = 8.0 * mcts_g.comm_costs[parent_id][task_id] / (mcts_g.bandwidth[parent_proc][processor] * log2(1 + global.P_BS * global.g_bs / double(global.Noise))) + 8.0 * mcts_g.comm_costs[parent_id][task_id] / MBS_trans_rate;
        else
            comm_time = 8.0 * mcts_g.comm_costs[parent_id][task_id] / MBS_trans_rate;

    }
    else {
        comm_time = 0;
    }
    return comm_time;
}

//------------------------------------------------------------------------//

void print_queue(Queue queue) {
    Element p;
    printf("Queue length:%d\n", queue->n);
    p = queue->head;
    printf("ready queue tasks id:");
    while (p) {
        printf("%d ", p->task->tag);
        p = p->next;
    }
    printf("\n");
}


Task delete_n_queue(Queue queue, int n) {
    int i = 0;
    Element p, prev = NULL;
    Task task;
    p = queue->head;
    if (n == 0) {
        queue->head = p->next;
        queue->n--;
        task = p->task;
        free(p);
        return task;
    }
    while (i < n) {
        prev = p;
        p = p->next;
        i++;
    }
    prev->next = p->next;
    queue->n--;
    task = p->task; free(p);
    return task;
}



Task delete_task_queue(Queue queue, int task_id) {
    Element p, prev;
    Task task;
    p = queue->head;
    if (p->task->tag == task_id) {
        queue->head = p->next;
        queue->n--;
        task = p->task;
        free(p);
        return task;
    }
    while (p) {
        prev = p;
        p = p->next;
        if (p && p->task->tag == task_id) {
            prev->next = p->next;
            queue->n--;
            task = p->task;
            free(p);
            return task;
        }
    }
    return NULL;
}

void insert_queue(Queue queue, Task task) {
    Element element;
    element = (Element)calloc(1, sizeof(struct _Element));
    element->task = task;
    element->next = queue->head;
    queue->head = element;
    queue->n++;
}

//The queue header adds the nb_tasks ready task
void update_queue(Queue queue, Task* task_ready, int nb_tasks) {
    int i;
    for (i = 0; i < nb_tasks; i++) {
        insert_queue(queue, task_ready[i]);
    }
}
//Initialize the queue. At this point, the queue does not contain any sub-tasks.
Queue init_queue() {
    Queue queue;
    queue = (Queue)calloc(1, sizeof(struct _Queue));
    queue->n = 0;
    queue->head = NULL;
    return queue;
}

void print_mcts_info(Tree tree, Queue queue) {
    Element p;
    p = queue->head;
    while (p)
    {
        printf("queue element %d\n", p->task->tag);
        p = p->next;
    }
    printf("tree %d:\n", tree->root->processor);
}

//Initializes child nodes based on the number of processors and queues
void init_children(Node node, Queue queue) {
    int i, j, k;
    node->nb_children = queue->n * global.n_processors;
    if (node->nb_children == 0) {
        node->children = NULL;
    }
    else {
        node->children = (Node*)calloc(node->nb_children, sizeof(Node));
        j = 0;
        for (i = 0; i < global.n_processors; i++) {
            Element p;
            p = queue->head;
            while (p) {
                node->children[j] = (Node)calloc(1, sizeof(struct _Node));
                node->children[j]->node_id = mcts_g.node_id++;
                node->children[j]->task_id = p->task->tag;
                node->children[j]->task = p->task;
                node->children[j]->processor = i;
                node->children[j]->Q = 0;
                node->children[j]->visits = 0;
                node->children[j]->nb_children = 0;
                node->children[j]->children = NULL;
                p = p->next;
                j++;
            }
        }
    }
}

//Initializes the root and its child nodes of the tree according to the dag
Task init_root_task(DAG dag)
{
    Task task;
    int i;
    task = (Task)calloc(1, sizeof(struct _Task));
    task->tag = 0;
    task->cost = 0;
    task->data_size = 0;
    task->priority = MAX_priority;
    task->nb_children = dag->nb_tasks_per_level[0];
    task->nb_parents = 0;
    task->nb_parents_r = 0;
    task->children = (Task*)calloc(task->nb_children, sizeof(Task));
    task->comm_costs = (double*)calloc(task->nb_children, sizeof(double));
    for (i = 0; i < dag->nb_tasks_per_level[0]; i++)
    {
        task->children[i] = dag->levels[0][i];
        task->comm_costs[i] = 0;
    }
    return task;
}


// Initialize the root node of the tree based on the DAG.
// The first level tasks of the DAG(nb_tasks_per_level[0]) are set as the child nodes of the root node.
// The child nodes of the root node are then added to the queue as ready nodes.
Tree init_tree(DAG dag, Queue queue) {
    int i;
    Tree tree;
    Task* task_ready;
    Node rootNode;
    //FILE* fp_process;
    //fp_process = fopen("makes.txt", "w");
    //fclose(fp_process);
    //printf("init_tree:1\n");
    //

    tree = (Tree)calloc(1, sizeof(struct _Tree));
    tree->depth = 1;
    tree->root = (Node)calloc(1, sizeof(struct _Node));
    rootNode = tree->root;
    rootNode->node_id = mcts_g.node_id++;
    rootNode->task_id = 0;
    rootNode->task = init_root_task(dag);
    rootNode->processor = 0;
    rootNode->Q = 0;
    rootNode->visits = 0;

    task_ready = (Task*)calloc(dag->nb_tasks_per_level[0], sizeof(Task));
    for (i = 0; i < dag->nb_tasks_per_level[0]; i++)
    {
        task_ready[i] = dag->levels[0][i];
    }
    update_queue(queue, task_ready, dag->nb_tasks_per_level[0]);// Add dag->nb_tasks_per_level[0] ready tasks to the front of the queue, i.e., task_ready
    init_children(tree->root, queue);
    mcts_g.nodetrace[mcts_g.nb_node] = rootNode;
    mcts_g.nb_node++;
    return tree;
}

//Select the child node with the largest uct value under this node
Node best_child(Node node) {
    int i;
    double uct, max, rp;
    Node bestNode = NULL;
    max = -9999999;
    rp = rand() / (RAND_MAX + 1.0);
    if (rp > 0) {
        for (i = 0; i < node->nb_children; i++) {
            int visits = node->children[i]->visits;
            if (visits != 0) {
                uct = -node->children[i]->Q / visits
                    - global.cp * sqrt(2 * log(node->visits) / visits);//UCB function
            }
            else {
                uct = INFINITY;
            }
            if (uct > max) {
                bestNode = node->children[i];
                max = uct;
            }
        }
    }
    else {
        i = rand() % node->nb_children;
        bestNode = node->children[i];
    }
    return bestNode;
}

//Treat the child nodes of the task as ready nodes and add these ready nodes to the queue
void update_ready_task(DAG dag, Task task, Queue queue)
{
    int i, nb_ready = 0;
    Task* task_ready = NULL;

    for (i = 0; i < task->nb_children; i++) // Traverse the child node of the task node
    {
        task->children[i]->nb_parents_r--; // The task has already been executed, so the precursor node of the child is minus one
        if (task->children[i]->nb_parents_r == 0) // If all of the child's precursor nodes are executed, the child is ready and enters the queue
        {
            Task* temp = (Task*)realloc(task_ready, (nb_ready + 1) * sizeof(Task));
            if (temp == NULL) {
                fprintf(stderr, "Error: Memory allocation failed for task_ready.\n");
                free(task_ready);
                return; // Exit if memory allocation fails
            }
            task_ready = temp;
            task_ready[nb_ready] = task->children[i];
            nb_ready++;
        }
    }

    if (nb_ready != 0)
    {
        update_queue(queue, task_ready, nb_ready); // A ready task (task_ready) is added to the queue header
    }

    // Free the allocated memory for task_ready after updating the queue
    free(task_ready);
}


double get_start_time(Task task, int processor) {
    int task_id, task_priority, i, j, parent_id, parent_proc;
    double comm_time;
    double start_time = 0;
    double finish_time = 0;
    task_id = task->tag;
    task_priority = task->priority;

    for (i = 0; i < task->nb_parents; i++)
    {

        parent_id = task->parents[i]->tag;
        parent_proc = mcts_g.scheduled[parent_id]->processor;

        comm_time = get_comm_time(parent_id, parent_proc, task_id, processor);
        start_time = start_time > (mcts_g.scheduled[parent_id]->EFT + comm_time) ? start_time : (mcts_g.scheduled[parent_id]->EFT + comm_time);
    }
    finish_time = start_time + task->comp_costs[processor];

    for (int k = 0; k < global.n + 1; k++)
    {
        i = mcts_g.scheduled_sort_on_per_server[processor][k];
        if (mcts_g.scheduled[i]->processor == processor)
        {

            if (finish_time > mcts_g.scheduled[i]->EST && start_time < mcts_g.scheduled[i]->EFT)//Jump-queue
            {
                if (mcts_g.scheduled[i]->priority > task_priority)
                {
                    start_time = mcts_g.scheduled[i]->EFT;
                    finish_time = start_time + task->comp_costs[processor];
                    // mcts_g.scheduled[task_id]->EST = start_time;
                    // mcts_g.scheduled[task_id]->EFT = finish_time;
                     //refresh_time_on_processor(task_id, processor);
                }
                else
                {
                    double tmp = mcts_g.scheduled[i]->EFT - mcts_g.scheduled[i]->EST;
                    double add_time = finish_time - mcts_g.scheduled[i]->EST;
                    mcts_g.scheduled[i]->EST = finish_time;
                    mcts_g.scheduled[i]->EFT = mcts_g.scheduled[i]->EST + tmp;
                    refresh_time_on_processor(i, mcts_g.scheduled[i]->processor);
                    refresh_time_on_succedtask(i, add_time);

                }

            }

        }
    }

    return start_time;
}
void refresh_time_on_processor(int task_id, int processor)
{
    for (int k = 0; k < global.n + 1; k++)
    {
        int j = mcts_g.scheduled_sort_on_per_server[processor][k];
        if (mcts_g.scheduled[j]->processor == processor && j != task_id)
        {
            if (mcts_g.scheduled[j]->EFT > mcts_g.scheduled[task_id]->EST && mcts_g.scheduled[j]->EST < mcts_g.scheduled[task_id]->EFT)
            {
                double tmp = mcts_g.scheduled[j]->EFT - mcts_g.scheduled[j]->EST;
                double add_time = mcts_g.scheduled[task_id]->EFT - mcts_g.scheduled[j]->EST;
                mcts_g.scheduled[j]->EST = mcts_g.scheduled[task_id]->EFT;
                mcts_g.scheduled[j]->EFT = mcts_g.scheduled[j]->EST + tmp;
                refresh_time_on_processor(j, mcts_g.scheduled[j]->processor);
                refresh_time_on_succedtask(j, add_time);

            }
        }

    }
}
void refresh_time_on_succedtask(int task_id, double add_time)
{
    Task task;
    task = mcts_g.scheduled[task_id]->task;
    for (int i = 0; i < task->nb_children; i++)
    {
        Task child = task->children[i];
        double comm_time;
        int task_processor = mcts_g.scheduled[task_id]->processor;
        int child_processor = mcts_g.scheduled[child->tag]->processor;
        comm_time = get_comm_time(task_id, task_processor, child->tag, child_processor);
        if (mcts_g.scheduled[child->tag]->EFT != 0
            && mcts_g.scheduled[child->tag]->EFT - mcts_g.scheduled[task_id]->EFT < comm_time)
        {
            double tmp = mcts_g.scheduled[child->tag]->EFT - mcts_g.scheduled[child->tag]->EST;
            mcts_g.scheduled[child->tag]->EST = mcts_g.scheduled[task_id]->EFT + comm_time;
            mcts_g.scheduled[child->tag]->EFT = mcts_g.scheduled[child->tag]->EST + tmp;
            refresh_time_on_processor(child->tag, mcts_g.scheduled[child->tag]->processor);
            refresh_time_on_succedtask(child->tag, add_time);

        }

    }
}


double find_start_time(Task task, int processor) {
    int task_id, task_priority, i, j, parent_id, parent_proc;
    double comm_time;
    double start_time = 0;
    double finish_time = 0;
    task_id = task->tag;
    task_priority = task->priority;

    for (i = 0; i < task->nb_parents; i++)
    {
        parent_id = task->parents[i]->tag;
        parent_proc = mcts_g.scheduled[parent_id]->processor;

        comm_time = get_comm_time(parent_id, parent_proc, task_id, processor);
        start_time = start_time > mcts_g.scheduled[parent_id]->EFT + comm_time ? start_time : mcts_g.scheduled[parent_id]->EFT + comm_time;
    }
    finish_time = start_time + task->comp_costs[processor]; // 
    for (i = 0; i < global.n + 1; i++)
    {
        if (mcts_g.scheduled[i]->processor == processor)
        {
            if (finish_time > mcts_g.scheduled[i]->EST && start_time < mcts_g.scheduled[i]->EFT)//Jump-queue
            {
                if (mcts_g.scheduled[i]->priority >= task_priority)
                {
                    start_time = mcts_g.scheduled[i]->EFT;
                    finish_time = start_time + task->comp_costs[processor];
                }
            }
        }
    }

    return start_time;
}

void update_scheduled_sort_on_processor(int task_id)
{

    int processor = mcts_g.scheduled[task_id]->processor;
    double est = mcts_g.scheduled[task_id]->EST;

    for (int i = global.n; i >= 0; i--)
    {
        int tmp_id = mcts_g.scheduled_sort_on_per_server[processor][i];
        if (est >= mcts_g.scheduled[tmp_id]->EST)
        {

            for (int j = 0; j <= i - 1; j++)
                mcts_g.scheduled_sort_on_per_server[processor][j] = mcts_g.scheduled_sort_on_per_server[processor][j + 1];
            mcts_g.scheduled_sort_on_per_server[processor][i] = task_id;
            break;
        }
    }


}


void update_scheduled_task(Task task, int processor) {
    double start_time;
    start_time = get_start_time(task, processor);
    mcts_g.scheduled[task->tag]->task = task;
    mcts_g.scheduled[task->tag]->priority = task->priority;
    mcts_g.scheduled[task->tag]->processor = processor;
    mcts_g.scheduled[task->tag]->EST = start_time;
    mcts_g.scheduled[task->tag]->EFT = start_time + task->comp_costs[processor];
    update_scheduled_sort_on_processor(task->tag);

}


void expand(Node node, Queue queue)
{
    init_children(node, queue);
}



void backup() {
    int i;
    double makespan;
    makespan = get_makespan();
    for (i = 0; i < mcts_g.nb_node; i++) {
        // printf("%lf\n", makespan);
        mcts_g.nodetrace[i]->Q += makespan;
        mcts_g.nodetrace[i]->visits++;
    }
}

void output_subtree(Node node, FILE* fp) {
    int p_nid, c_nid, p_tid, c_tid;
    double uct;
    if (node != NULL) {
        fprintf(fp, "%d [label = \"Q=%.2lf\nN=%d\n Q/N=%.2lf\"]\n",
            node->node_id, node->Q, node->visits, node->Q / node->visits);
        p_nid = node->node_id; p_tid = node->task_id;
        for (int i = 0; i < node->nb_children; i++) {
            c_nid = node->children[i]->node_id; c_tid = node->children[i]->task_id;
            uct = -node->children[i]->Q / node->children[i]->visits
                - CP * sqrt(2 * log(node->visits) / node->children[i]->visits);
            fprintf(fp, "%d -> %d ", p_nid, c_nid);
            fprintf(fp, "[label = \"t=%d, p=%d, uct=%.2lf\"]\n", c_tid - 1, node->children[i]->processor, uct);
            output_subtree(node->children[i], fp);
        }
    }
}


void output_tree(Tree tree) {
    FILE* fp;
    Node node;
    int p_id, c_id;
    errno_t err = fopen_s(&fp, "tree.dot", "w");
    //fp = fopen("tree.dot", "w");
    fprintf(fp, "digraph G {");
    node = tree->root;
    if (node != NULL) {
        output_subtree(node, fp);
    }
    fprintf(fp, "}\n");
    fclose(fp);
}

// Search the tree
Node tree_policy(DAG dag, Tree tree, Queue queue) {
    Node node, child;
    Task task;

    //FILE* fp_process;
    //errno_t err = fopen_s(&fp_process, "makes.txt", "a");
    //fp_process = fopen("makes.txt", "a");

    node = tree->root;
    while (node->children != NULL) {// Traverses all nodes

        child = best_child(node);// Select the optimal child node for the current node

        task = delete_task_queue(queue, child->task_id);// Select the optimal child node for the current node

        update_scheduled_task(task, child->processor);// Update the task status, recording the start time and completion time

        update_ready_task(dag, task, queue);// Update the ready node and add it to the queue

        mcts_g.nodetrace[mcts_g.nb_node] = child;// Record the child in the path that has been traveled
        mcts_g.nb_node++;
        node = child;

    }

    //fclose(fp_process);
    //printf("tree_policy:8\n");
    return node;
}


//Get the processor that minimizes the time it takes for the task to finish
int get_best_processor(Task task) {
    int i, bestp = 0;
    double est, eft, min = INFINITY;
    for (i = 0; i < global.n_processors; i++) {
        // printf("get_best_processer\n");
        est = find_start_time(task, i);
        eft = est + task->comp_costs[i];
        if (min > eft) {
            min = eft;
            bestp = i;
        }
    }
    return bestp;
}

//Determine the assignment of the leafnode
void default_random_policy(DAG dag, Node leafnode, Queue queue) {
    Task task, child;
    int processor;
    double rp;

    //FILE* fp_process;
    //errno_t err = fopen_s(&fp_process, "makes.txt", "a");
    //fp_process = fopen("makes.txt", "a");

    //    printf("default_random_policy:1\n");
    task = leafnode->task;
    while (queue->n != 0)
    {
        child = delete_n_queue(queue, rand() % queue->n);
        update_scheduled_task(child, get_best_processor(child));
        //rp = rand() / (RAND_MAX + 1.0);

        //if (rp > 0.2)
        //{
        //   update_scheduled_task(child, get_best_processor(child), fp_process);
        //}
        //else 
        //{
        //    update_scheduled_task(child, rand() % global.n_processors, fp_process);
        //}
        //update_scheduled_task(child, get_best_processor(child));

        update_ready_task(dag, child, queue);
    }

    //fclose(fp_process);
}

// Obtain the finish time of each task
double get_makespan() {
    int i;
    double eft, max = -1;
    for (i = 0; i < global.n + 1; i++) {
        eft = mcts_g.scheduled[i]->EFT;
        if (eft > max) {
            max = eft;
        }
    }
    return max;
}

double get_avg_satisfaction(DAG dag)
{
    double sum_satisfaction = 0;
    for (int i = 0; i < global.nb_users; i++)
    {
        Task task = dag->levels[dag->nb_levels - 1][i];
        if (global.user_satisfaction[i] <= 3)
            if (mcts_g.scheduled[task->tag]->EFT <= global.user_satisfaction[i])
                sum_satisfaction += log(1 + global.user_satisfaction[i] - mcts_g.scheduled[task->tag]->EFT) / 2.0;
            else
                sum_satisfaction += -1;
        else if (global.user_satisfaction[i] > 3 && global.user_satisfaction[i] <= 6)
            if (mcts_g.scheduled[task->tag]->EFT <= global.user_satisfaction[i])
                sum_satisfaction += 1;
            else
                sum_satisfaction += 1 * exp(-(mcts_g.scheduled[task->tag]->EFT - global.user_satisfaction[i]));
        else
            sum_satisfaction += 1;
    }
    return sum_satisfaction / global.nb_users;
}


double get_speedup(DAG dag) {

    int best_performance_processor = 0;
    for (int i = 1; i < global.n_processors; i++)
        if (mcts_g.processor_performance[i] > mcts_g.processor_performance[best_performance_processor])
            best_performance_processor = i;
    for (int i = 0; i < dag->nb_levels; i++)
        for (int j = 0; j < dag->nb_tasks_per_level[i]; j++)
        {
            update_scheduled_task(dag->levels[i][j], best_performance_processor);
        }
    double makespan = get_makespan();
    mcts_g.scheduled = (Scheduled*)calloc(global.n + 1, sizeof(Scheduled));
    for (int i = 0; i < global.n + 1; i++)
    {
        mcts_g.scheduled[i] = (Scheduled)calloc(1, sizeof(struct _scheduled));
    }
    mcts_g.scheduled_sort_on_per_server = (int**)calloc(global.n_processors, sizeof(int*));
    for (int i = 0; i < global.n_processors; i++)
        mcts_g.scheduled_sort_on_per_server[i] = (int*)calloc(global.n + 1, sizeof(int));

    return makespan;
}

//Maximize reward
void UCTsearch(DAG dag, Tree tree, Queue queue) {
    int i, j, k, nb_simulation = 0;
    Node leaf_node;
    Task* task_ready, task;
    double temp_satisfaction, max_satisfaction = -1000;
    //FILE* fp_ms;
    //errno_t err = fopen_s(&fp_ms, "makespan.txt", "a");
    //fp_ms = fopen("makespan.txt", "w");

    //FILE* fp_process;
    //errno_t err_fp = fopen_s(&fp_process, "makes.txt", "a");
    //fp_process = fopen("makes.txt", "a");

    //fprintf(fp_ms, "%d %d %d %d \n", global.n, global.n_processors, global.n_simulations, global.cp);// These values are already defined in 'daggen_commons.h'
    int* best_decision = (int*)calloc(global.n + 1, sizeof(int));
    double* best_decision_eft = (double*)calloc(global.nb_users, sizeof(double));

    while (nb_simulation < global.n_simulations)
    {
        leaf_node = tree_policy(dag, tree, queue);//Search the tree from the root node, until reach to the leafnode
        if (leaf_node->visits == 0)
        {
            //If this leafnode has not been visited
            default_random_policy(dag, leaf_node, queue);//Simulate the remaining nodes in the queue
        }
        else
        {
            //If this leafnode has been visited
            expand(leaf_node, queue);
            if (leaf_node->nb_children != 0)
            {
                task = delete_task_queue(queue, leaf_node->children[0]->task_id);
                update_scheduled_task(task, leaf_node->children[0]->processor);
                update_ready_task(dag, task, queue);
                mcts_g.nodetrace[mcts_g.nb_node] = leaf_node->children[0];
                mcts_g.nb_node++;
                default_random_policy(dag, leaf_node->children[0], queue);
            }
            else //leaf_node->nb_children = 0
            {
                default_random_policy(dag, leaf_node, queue);
            }
        }
        backup();
        temp_satisfaction = get_avg_satisfaction(dag);
        //fprintf(fp_ms, "%lf %lf\n", temp_satisfaction, mcts_g.nodetrace[0]->Q / mcts_g.nodetrace[0]->visits);
        if (max_satisfaction < temp_satisfaction)// Update the current best strategy 
        {
            max_satisfaction = temp_satisfaction;
            result.mcts_makespan = get_makespan();
            //for (int i = 1; i < global.n + 1; i++)
            //{
            //    best_decision[mcts_g.scheduled[i]->task->tag] = mcts_g.scheduled[i]->processor;
            //}
            for (int i = 0; i < global.nb_users; i++)
            {
                Task task = dag->levels[dag->nb_levels - 1][i];
                best_decision_eft[i] = mcts_g.scheduled[task->tag]->EFT;
            }
            //for (int i = 1; i < global.n + 1; i++)
            //{
            //printf("id, processor, priority, EST,EFT: %d %d, %d, %lf, %lf\n", mcts_g.scheduled[i]->task->tag, mcts_g.scheduled[i]->processor, mcts_g.scheduled[i]->priority, mcts_g.scheduled[i]->EST, mcts_g.scheduled[i]->EFT);
            //}
            //printf("\n");
            //for(i=0;i<mcts_g.nb_node;i++) 
           // printf("%d ",mcts_g.nodetrace[i]->task_id); 
            //printf("\n");
        }


        k = 0;
        for (i = 1; i < dag->nb_levels; i++)
            for (j = 0; j < dag->nb_tasks_per_level[i]; j++) {
                dag->levels[i][j]->nb_parents_r = global.nb_parents[k];
                // printf("i, j, k, global.nb_parents[k]: %d %d %d %d \n", i, j, k, global.nb_parents[k]);
                k++;
            }

        //Reset for the next iteration
        mcts_g.nb_node = 1;

        //Initialize the  queue
        task_ready = (Task*)calloc(dag->nb_tasks_per_level[0], sizeof(Task));
        for (i = 0; i < dag->nb_tasks_per_level[0]; i++) {
            task_ready[i] = dag->levels[0][i];
        }
        update_queue(queue, task_ready, dag->nb_tasks_per_level[0]);

        // Allocate and initialize mcts_g.scheduled
        for (i = 0; i < global.n + 1; i++) {
            if (mcts_g.scheduled[i] != NULL) {
                free(mcts_g.scheduled[i]);
            }
        }
        free(mcts_g.scheduled);

        mcts_g.scheduled = (Scheduled*)calloc(global.n + 1, sizeof(Scheduled));
        for (i = 0; i < global.n + 1; i++) {
            mcts_g.scheduled[i] = (Scheduled)calloc(1, sizeof(struct _scheduled));
        }

        // Allocate and initialize mcts_g.scheduled_sort_on_per_server
        for (i = 0; i < global.n_processors; i++) {
            if (mcts_g.scheduled_sort_on_per_server[i] != NULL) {
                free(mcts_g.scheduled_sort_on_per_server[i]);
            }
        }
        free(mcts_g.scheduled_sort_on_per_server);

        mcts_g.scheduled_sort_on_per_server = (int**)calloc(global.n_processors, sizeof(int*));
        for (i = 0; i < global.n_processors; i++) {
            mcts_g.scheduled_sort_on_per_server[i] = (int*)calloc(global.n + 1, sizeof(int));
        }

        nb_simulation++;

    }
    //output_tree(tree);
    //fclose(fp_process);

    //fclose(fp_ms);
    result.mcts_satisfaction = max_satisfaction;
    result.mcts_speedup = get_speedup(dag) / result.mcts_makespan;
    free(queue);
    free(tree->root);
    free(tree);
    //free(task_ready);
    free(best_decision);
    free(best_decision_eft);

}

