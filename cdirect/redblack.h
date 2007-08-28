#ifndef REDBLACK_H
#define REDBLACK_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef enum { RED, BLACK } rb_color;
typedef struct rb_node_s {
     struct rb_node_s *p, *r, *l; /* parent, right, left */
     int k; /* key/data ... for DIRECT, an index into our hyperrect array */
     rb_color c;
} rb_node;

typedef int (*rb_compare)(int k1, int k2, void *c_data);

typedef struct {
     rb_compare compare; void *c_data;
     rb_node *root;
     int N; /* number of nodes */

     /* in our application, we can optimize memory allocation because
	we never delete two nodes in a row (we always add a node after
	deleting)... or rather, we never delete but the value of
        the key sometimes changes.  ... this means we can just
	allocate a linear, exponentially growing stack (nodes) of
	nodes, and don't have to worry about holes in the stack */
     rb_node *nodes; /* allocated data of nodes, in some order */
     int Nalloc; /* number of allocated nodes */
     rb_node nil; /* explicit node for NULL nodes, for convenience */
} rb_tree;

extern int rb_tree_init(rb_tree *t, rb_compare compare, void *c_data);
extern void rb_tree_destroy(rb_tree *t);
extern int rb_tree_insert(rb_tree *t, int k);
extern int rb_tree_check(rb_tree *t);
extern rb_node *rb_tree_find(rb_tree *t, int k);
extern rb_node *rb_tree_find_exact(rb_tree *t, int k);
extern rb_node *rb_tree_resort(rb_tree *t, rb_node *n);
extern rb_node *rb_tree_min(rb_tree *t);
extern rb_node *rb_tree_max(rb_tree *t);
extern rb_node *rb_tree_succ(rb_tree *t, rb_node *n);
extern rb_node *rb_tree_pred(rb_tree *t, rb_node *n);

/* To change a key, use rb_tree_find+resort.  Removing a node
   currently wastes memory unless you change the allocation scheme
   in redblack.c */
extern rb_node *rb_tree_remove(rb_tree *t, rb_node *n);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
