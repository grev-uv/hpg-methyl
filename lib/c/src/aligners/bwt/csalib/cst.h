/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
//////////////////////////////////////////////////////////////////////////////////////
//   functions for suffix trees
//////////////////////////////////////////////////////////////////////////////////////
#ifndef _CST_H_
#define _CST_H_
#include "csa.h"

typedef struct {
  CSA *csa;
  i64 depth;  // depth of the node
  i64 depth2;
  i64 l,r;    // range [l,r] of the string
} cst_node;

cst_node cst_root(CSA *csa);
int cst_isnull(cst_node node);
int cst_isleaf(cst_node node);
int cst_eq(cst_node n1, cst_node n2);
int cst_isunary(cst_node node);
uchar *cst_pathlabel(cst_node node);
i64 cst_depth(cst_node node);
cst_node cst_child(cst_node node, int c);
cst_node cst_firstchild(cst_node node);
cst_node cst_nextchild(cst_node node, cst_node child);
cst_node cst_weiner_link(cst_node node, int c);
cst_node cst_canonize(cst_node node);
cst_node cst_parent(cst_node node);
cst_node cst_suflink(cst_node node);


#endif // _CST_H_
