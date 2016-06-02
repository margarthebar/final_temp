#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct vary_node{
  char name[128];
  double value;
  struct vary_node *next;
};

char name1[128] = "first";
char name2[128] = "second";

struct vary_node ** second_pass(){
  struct vary_node *knobs[10];
  int i;
  for(i=0; i<10; i++){
    struct vary_node *next_node = malloc (sizeof(struct vary_node));
    sprintf(next_node->name,"%s",name1);
    next_node->value = 42;
    next_node->next = NULL;
    knobs[i] = next_node;
  }
  int j;
  for(j=0; j<10; i++){
    struct vary_node *next_node = malloc (sizeof(struct vary_node));
    sprintf(next_node->name,"%s",name2);
    next_node->value = 42;
    next_node->next = NULL;
    knobs[i]->next = next_node;
  }
  return knobs;
}

void main(){
  struct vary_node *knobs[10];
  knobs = second_pass();
  printf("done!\n");
}
