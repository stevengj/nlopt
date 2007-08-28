/* simple implementation of red-black trees optimized for use with DIRECT */

#include <stddef.h>
#include <stdlib.h>
#include "redblack.h"

int rb_tree_init(rb_tree *t, rb_compare compare, void *c_data) {
     t->compare = compare; t->c_data = c_data;
     t->nil.c = BLACK; t->nil.l = t->nil.r = t->nil.p = &t->nil; t->nil.k = -1;
     t->root = &t->nil;
     t->N = 0;
     t->Nalloc = 100; /* allocate some space to start with */
     t->nodes = (rb_node*) malloc(sizeof(rb_node) * t->Nalloc);
     return t->nodes != NULL;
}

void rb_tree_destroy(rb_tree *t)
{
     t->root = 0; t->N = 0; t->Nalloc = 0;
     free(t->nodes); t->nodes = 0;
}

/* in our application, we can optimize memory allocation because
   we never delete two nodes in a row (we always add a node after
   deleting)... or rather, we never delete but the value of
   the key sometimes changes.  ... this means we can just
   allocate a linear, exponentially growing stack (nodes) of
   nodes, and don't have to worry about holes in the stack ...
   otherwise, alloc1 should be replaced by an implementation that
   malloc's each node separately */
static rb_node *alloc1(rb_tree *t, int k)
{
     rb_node *nil = &t->nil;
     rb_node *n;
     if (t->Nalloc == t->N) { /* grow allocation */
	  rb_node *old_nodes = t->nodes;
	  ptrdiff_t change;
	  int i;
	  t->Nalloc = 2*t->Nalloc + 1;
	  t->nodes = (rb_node*) realloc(t->nodes, sizeof(rb_node) * t->Nalloc);
	  if (!t->nodes) return NULL;
	  change = t->nodes - old_nodes;
	  if (t->root != nil) t->root += change;
	  for (i = 0; i < t->N; ++i) { /* shift all pointers, ugh */
	       if (t->nodes[i].p != nil) t->nodes[i].p += change;
	       if (t->nodes[i].r != nil) t->nodes[i].r += change;
	       if (t->nodes[i].l != nil) t->nodes[i].l += change;
	  }
     }
     n = t->nodes + t->N++;
     n->k = k;
     n->p = n->l = n->r = nil;
     return n;
}

static void rotate_left(rb_node *p, rb_tree *t)
{
     rb_node *nil = &t->nil;
     rb_node *n = p->r; /* must be non-NULL */
     p->r = n->l;
     n->l = p;
     if (p->p != nil) {
	  if (p == p->p->l) p->p->l = n;
	  else p->p->r = n;
     }
     else
	  t->root = n;
     n->p = p->p;
     p->p = n;
     if (p->r != nil) p->r->p = p;
}

static void rotate_right(rb_node *p, rb_tree *t)
{
     rb_node *nil = &t->nil;
     rb_node *n = p->l; /* must be non-NULL */
     p->l = n->r;
     n->r = p;
     if (p->p != nil) {
	  if (p == p->p->l) p->p->l = n;
	  else p->p->r = n;
     }
     else
	  t->root = n;
     n->p = p->p;
     p->p = n;
     if (p->l != nil) p->l->p = p;
}

static void insert_node(rb_tree *t, rb_node *n)
{
     rb_node *nil = &t->nil;
     rb_compare compare = t->compare;
     void *c_data = t->c_data;
     int k = n->k;
     rb_node *p = t->root;
     n->c = RED;
     if (p == nil) {
	  t->root = n;
	  n->c = BLACK;
	  return;
     }
     /* insert (RED) node into tree */
     while (1) {
	  if (compare(k, p->k, c_data) <= 0) { /* k <= p->k */
	       if (p->l != nil)
		    p = p->l;
	       else {
		    p->l = n;
		    n->p = p;
		    break;
	       }
	  }
	  else {
	       if (p->r != nil)
		    p = p->r;
	       else {
		    p->r = n;
		    n->p = p;
		    break;
	       }
	  }
     }
 fixtree:
     if (n->p->c == RED) { /* red cannot have red child */
	  rb_node *u = p == p->p->l ? p->p->r : p->p->l;
	  if (u != nil && u->c == RED) {
	       p->c = u->c = BLACK;
	       n = p->p;
	       if ((p = n->p) != nil) {
		    n->c = RED;
		    goto fixtree;
	       }
	  }
	  else {
	       if (n == p->r && p == p->p->l) {
		    rotate_left(p, t);
		    p = n; n = n->l;
	       }
	       else if (n == p->l && p == p->p->r) {
		    rotate_right(p, t);
		    p = n; n = n->r;
	       }
	       p->c = BLACK;
	       p->p->c = RED;
	       if (n == p->l && p == p->p->l)
		    rotate_right(p->p, t);
	       else if (n == p->r && p == p->p->r)
		    rotate_left(p->p, t);
	  }
	      
     }
}

int rb_tree_insert(rb_tree *t, int k)
{
     rb_node *n = alloc1(t, k);
     if (!n) return 0;
     insert_node(t, n);
     return 1;
}

static int check_node(rb_node *n, int *nblack, rb_tree *t)
{
     rb_node *nil = &t->nil;
     int nbl, nbr;
     rb_compare compare = t->compare;
     void *c_data = t->c_data;
     if (n == nil) { *nblack = 0; return 1; }
     if (n->r != nil && n->r->p != n) return 0;
     if (n->r != nil && compare(n->r->k, n->k, c_data) < 0)
	  return 0;
     if (n->l != nil && n->l->p != n) return 0;
     if (n->l != nil && compare(n->l->k, n->k, c_data) > 0)
	  return 0;
     if (n->c == RED) {
	  if (n->r != nil && n->r->c == RED) return 0;
	  if (n->l != nil && n->l->c == RED) return 0;
     }
     if (!(check_node(n->r, &nbl, t) && check_node(n->l, &nbr, t))) 
	  return 0;
     if (nbl != nbr) return 0;
     *nblack = nbl + n->c == BLACK;
     return 1;
}
int rb_tree_check(rb_tree *t)
{
     rb_node *nil = &t->nil;
     int nblack;
     if (nil->c != BLACK) return 0;
     if (t->root == nil) return 1;
     if (t->root->c != BLACK) return 0;
     return check_node(t->root, &nblack, t);
}

rb_node *rb_tree_find(rb_tree *t, int k)
{
     rb_node *nil = &t->nil;
     rb_compare compare = t->compare;
     void *c_data = t->c_data;
     rb_node *p = t->root;
     while (p != nil) {
	  int comp = compare(k, p->k, c_data);
	  if (!comp) return p;
	  p = comp <= 0 ? p->l : p->r;
     }
     return NULL;
}

/* like rb_tree_find, but guarantees that returned node n will have
   n->k == k (may not be true above if compare(k,k') == 0 for some k != k') */
rb_node *rb_tree_find_exact(rb_tree *t, int k)
{
     rb_node *nil = &t->nil;
     rb_compare compare = t->compare;
     void *c_data = t->c_data;
     rb_node *p = t->root;
     while (p != nil) {
	  int comp = compare(k, p->k, c_data);
	  if (!comp) break;
	  p = comp <= 0 ? p->l : p->r;
     }
     if (p == nil)
	  return NULL;
     while (p->l != nil && !compare(k, p->l->k, c_data)) p = p->l;
     if (p->l != nil) p = p->l;
     do {
	  if (p->k == k) return p;
	  p = rb_tree_succ(t, p);
     } while (p && compare(p->k, k, c_data) <= 0);
     return NULL;
}

/* find greatest point in subtree p that is <= k */
static rb_node *find_le(rb_node *p, int k, rb_tree *t)
{
     rb_node *nil = &t->nil;
     rb_compare compare = t->compare;
     void *c_data = t->c_data;
     while (p != nil) {
	  if (compare(p->k, k, c_data) <= 0) { /* p->k <= k */
	       rb_node *r = find_le(p->r, k, t);
	       if (r) return r;
	       else return p;
	  }
	  else /* p->k > k */
	       p = p->l;
     }
     return NULL; /* k < everything in subtree */
}

/* find greatest point in t <= k */
rb_node *rb_tree_find_le(rb_tree *t, int k)
{
     return find_le(t->root, k, t);
}

/* find least point in subtree p that is > k */
static rb_node *find_gt(rb_node *p, int k, rb_tree *t)
{
     rb_node *nil = &t->nil;
     rb_compare compare = t->compare;
     void *c_data = t->c_data;
     while (p != nil) {
	  if (compare(p->k, k, c_data) > 0) { /* p->k > k */
	       rb_node *l = find_gt(p->l, k, t);
	       if (l) return l;
	       else return p;
	  }
	  else /* p->k <= k */
	       p = p->r;
     }
     return NULL; /* k >= everything in subtree */
}

/* find least point in t > k */
rb_node *rb_tree_find_gt(rb_tree *t, int k)
{
     return find_gt(t->root, k, t);
}

rb_node *rb_tree_min(rb_tree *t)
{
     rb_node *nil = &t->nil;
     rb_node *n = t->root;
     while (n != nil && n->l != nil)
	  n = n->l;
     return(n == nil ? NULL : n);
}

rb_node *rb_tree_max(rb_tree *t)
{
     rb_node *nil = &t->nil;
     rb_node *n = t->root;
     while (n != nil && n->r != nil)
	  n = n->r;
     return(n == nil ? NULL : n);
}

rb_node *rb_tree_succ(rb_tree *t, rb_node *n)
{
     rb_node *nil = &t->nil;
     if (n->r == nil) {
	  rb_node *prev;
	  do {
	       prev = n;
	       n = n->p;
	  } while (prev == n->r && n != nil);
	  return n == nil ? NULL : n;
     }
     else {
	  n = n->r;
	  while (n->l != nil)
	       n = n->l;
	  return n;
     }
}

rb_node *rb_tree_pred(rb_tree *t, rb_node *n)
{
     rb_node *nil = &t->nil;
     if (n->l == nil) {
	  rb_node *prev;
	  do {
	       prev = n;
	       n = n->p;
	  } while (prev == n->l && n != nil);
	  return n == nil ? NULL : n;
     }
     else {
	  n = n->l;
	  while (n->r != nil)
	       n = n->r;
	  return n;
     }
}

rb_node *rb_tree_remove(rb_tree *t, rb_node *n)
{
     rb_node *nil = &t->nil;
     rb_node *m;
     if (n->l != nil && n->r != nil) {
	  rb_node *lmax = n->l;
	  while (lmax->r != nil) lmax = lmax->r;
	  n->k = lmax->k;
	  n = lmax;
     }
     m = n->l != nil? n->l : n->r;
     if (n->p != nil) {
	  if (n->p->r == n) n->p->r = m;
	  else n->p->l = m;
     }
     else
	  t->root = m;
     m->p = n->p;
     if (n->c == BLACK) {
	  if (m->c == RED)
	       m->c = BLACK;
	  else {
	  deleteblack:
	       if (m->p != nil) {
		    rb_node *s = m == m->p->l ? m->p->r : m->p->l;
		    if (s->c == RED) {
			 m->p->c = RED;
			 s->c = BLACK;
			 if (m == m->p->l) rotate_left(m->p, t);
			 else rotate_right(m->p, t);
			 s = m == m->p->l ? m->p->r : m->p->l;
		    }
		    if (m->p->c == BLACK && s->c == BLACK
			&& s->l->c == BLACK && s->r->c == BLACK) {
			 if (s != nil) s->c = RED;
			 m = m->p;
			 goto deleteblack;
		    }
		    else if (m->p->c == RED && s->c == BLACK &&
			     s->l->c == BLACK && s->r->c == BLACK) {
			 if (s != nil) s->c = RED;
			 m->p->c = BLACK;
		    }
		    else {
			 if (m == m->p->l && s->c == BLACK &&
			     s->l->c == RED && s->r->c == BLACK) {
			      s->c = RED;
			      s->l->c = BLACK;
			      rotate_right(s, t);
			      s = m == m->p->l ? m->p->r : m->p->l;
			 }
			 else if (m == m->p->r && s->c == BLACK &&
				  s->r->c == RED && s->l->c == BLACK) {
			      s->c = RED;
			      s->r->c = BLACK;
			      rotate_left(s, t);
			      s = m == m->p->l ? m->p->r : m->p->l;
			 }
			 s->c = m->p->c;
			 m->p->c = BLACK;
			 if (m == m->p->l) {
			      s->r->c = BLACK;
			      rotate_left(m->p, t);
			 }
			 else {
			      s->l->c = BLACK;
			      rotate_right(m->p, t);
			 }
		    }
	       }
	  }
     }
     return n; /* the node that was deleted may be different from initial n */
}

rb_node *rb_tree_resort(rb_tree *t, rb_node *n)
{
     int k = n->k;
     n = rb_tree_remove(t, n);
     n->p = n->l = n->r = &t->nil;
     n->k = k; /* n may have changed during remove */
     insert_node(t, n);
     return n;
}
