/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include "cst.h"

cst_node cst_root(CSA *csa)
{
  cst_node node;

  node.csa = csa;
  node.depth = 0;
  node.depth2 = 0;
  node.l = 0;
  node.r = csa->n;
  
  return node;
}

uchar *cst_pathlabel(cst_node node)
{
  uchar *buf;
  buf = (uchar *) mymalloc(node.depth+1);
  node.csa->substring(buf, node.csa, node.l, node.depth);
  buf[node.depth] = 0;
  return buf;
}

i64 cst_depth(cst_node node)
{
  return node.depth;
}

int cst_isnull(cst_node node)
{
  return (node.depth == -1);
}

int cst_isleaf(cst_node node)
{
  return (node.l == node.r);
}

int cst_eq(cst_node n1, cst_node n2)
{
  if (n1.csa != n2.csa) {
    return 0;
  }
  if (n1.depth != n2.depth) return 0;
  if (n1.l != n2.l) return 0;
  if (n1.r != n2.r) return 0;
  return 1;
}


int cst_isunary(cst_node node)
{
  int c1,c2;
  i64 l, r;
  i64 i, depth;
  CSA *csa;

  csa = node.csa;
  l = node.l;  r = node.r;
  depth = node.depth;

  for (i=0; i<depth; i++) {
    l = csa->psi(csa, l);
    r = csa->psi(csa, r);
  }

  c1 = csa->head(csa, l);
  c2 = csa->head(csa, r);

  return (c1 == c2);
}

cst_node cst_child(cst_node node, int c)
{
  cst_node child;
  i64 ll, rr;
  i64 len;
  
  child.csa = node.csa;

  ll = node.l;  rr = node.r;
  len = csa_search_r(node.depth, c, node.csa, &ll, &rr);

  if (len == node.depth+1) {
    child.depth = node.depth+1;
    child.depth2 = node.depth+1;
    child.l = ll;
    child.r = rr;
  } else {
    child.depth = -1;
    child.depth2 = -1;
    child.l = 1;
    child.r = 0;
  }
  return child;
}

cst_node cst_firstchild(cst_node node)
{
  int c;
  i64 l;
  i64 i, depth;
  CSA *csa;
  cst_node newnode;

  csa = node.csa;
  l = node.l;
  depth = node.depth;

  for (i=0; i<depth; i++) {
    l = csa->psi(csa, l);
  }

  c = csa->head(csa, l);
  if (c == -1) {
    newnode.csa = csa;
    newnode.depth = depth+1;
    newnode.depth2 = depth+1;
    newnode.l = newnode.r = node.l;
    return newnode;
  } else {
    return cst_child(node, csa->AtoC[c]);
  }
}

cst_node cst_nextchild(cst_node node, cst_node child)
{
  int c;
  i64 l;
  i64 i, depth;
  CSA *csa;
  cst_node newnode;

  csa = node.csa;
  l = child.r+1;
  depth = node.depth;

  newnode.csa = 0;
  if (l > node.r) {
    newnode.depth = -1;
    newnode.depth2 = -1;
    newnode.l = 1;
    newnode.r = 0;
    return newnode;
  }

  for (i=0; i<depth; i++) {
    l = csa->psi(csa, l);
  }

  c = csa->head(csa, l);
  if (c == -1) {
    newnode.csa = csa;
    newnode.depth = depth+1;
    newnode.depth2 = depth+1;
    newnode.l = newnode.r = node.l;
    return newnode;
  } else {
    return cst_child(node, csa->AtoC[c]);
  }
}

cst_node cst_weiner_link(cst_node node, int c)
{
  CSA *csa;
  cst_node newnode;
  i64 ll, rr;

  
  csa = newnode.csa = node.csa;

  ll = node.l;  rr = node.r;
  if (csa->searchsub(c, node.csa, &ll, &rr) != 0) {
    newnode.depth = -1;
    newnode.depth2 = -1;
    newnode.l = 1;
    newnode.r = 0;
  } else {
    newnode.depth = node.depth+1;
    newnode.depth2 = node.depth+1;
    newnode.l = ll;
    newnode.r = rr;
  }
  return newnode;

}

cst_node cst_canonize(cst_node node)
{
  int c1,c2;
  cst_node newnode;
  i64 l, r;
  i64 i, depth;
  CSA *csa;

  newnode = node;

  csa = node.csa;
  l = node.l;  r = node.r;
  depth = node.depth;

  for (i=0; i<depth; i++) {
    l = csa->psi(csa, l);
    r = csa->psi(csa, r);
  }

  while (1) {
    c1 = csa->head(csa, l);
    c2 = csa->head(csa, r);
    if (c1 != c2) break;
    l = csa->psi(csa, l);
    r = csa->psi(csa, r);
    depth++;
  }
  
  newnode.depth = depth;
  return newnode;
}

cst_node cst_parent(cst_node node)
{
  CSA *csa;
  uchar *label;
  cst_node parent;
  i64 l, r;


  csa = node.csa;
  parent.csa = csa;
  parent.depth = node.depth-1;
  parent.depth2 = node.depth-1;

  label = cst_pathlabel(node);
  l = node.l;  r = node.r;
  if (csa->search(label, parent.depth, csa, &l, &r) != parent.depth) {
    printf("cst_parent: ???\n");
  }
  parent.l = l;  parent.r = r;

  free(label);
  return parent;
}

cst_node cst_suflink(cst_node node)
{
  CSA *csa;
  i64 l, r;
  cst_node newnode;
  
  csa = newnode.csa = node.csa;

  l = node.l;  r = node.r;
  newnode.l = csa->psi(csa, l);
  newnode.r = csa->psi(csa, r);
  newnode.depth = node.depth;
  return cst_parent(newnode);
}
