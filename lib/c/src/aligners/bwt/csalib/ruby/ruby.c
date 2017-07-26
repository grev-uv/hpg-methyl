/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <ruby.h>
#include "../csa.h"

static CSA *
csa_ptr(VALUE obj)
{
    Check_Type(obj, T_DATA);
    return DATA_PTR(obj);
}

static void
csa_mark(void *ptr)
{
    return;
}

static void
csa_free(void *ptr)
{
    if (ptr) {
      xfree(ptr);
    }
}

static VALUE
rb_csa_alloc_func(VALUE klass)
{
    CSA *sa;
    Data_Make_Struct(klass, CSA, csa_mark, csa_free, sa);
}

static VALUE
rb_csa_initialize(VALUE self, VALUE idx, VALUE psi)
{
    CSA *sa = csa_ptr(self);
    char *fname[2];

    fname[0] = StringValueCStr(idx);
    fname[1] = StringValueCStr(psi);
    csa_read(sa, 2, fname);

    return self;
}

static VALUE
rb_csa_getn(VALUE self)
{
    CSA *sa = csa_ptr(self);
    return LONG2FIX(sa->n);
}

static VALUE
rb_csa_getm(VALUE self)
{
    CSA *sa = csa_ptr(self);
    return INT2FIX(sa->m);
}

static VALUE
rb_csa_charset(VALUE self)
{
    CSA *sa = csa_ptr(self);
    VALUE charset;
    int i,m;

    m = sa->m;
    charset = rb_ary_new();
    for (i=0; i<m; i++) {
      rb_ary_push(charset, INT2FIX(sa->AtoC[i]));
    }
    return charset;
}

static VALUE
rb_csa_atoc(VALUE self, VALUE oi)
{
    CSA *sa = csa_ptr(self);
    int i, j;

    i = FIX2INT(oi);
    if (i < 0 || i >= sa->m) {
      return Qnil;
    }
    j = sa->AtoC[i];

    return INT2FIX(j);
}

static VALUE
rb_csa_ctoa(VALUE self, VALUE oi)
{
    CSA *sa = csa_ptr(self);
    int i, j;

    i = FIX2INT(oi);
    if (i < 0 || i >= SIGMA) {
      return Qnil;
    }
    j = sa->CtoA[i];

    return INT2FIX(j);
}

static VALUE
rb_csa_bw(VALUE self, VALUE oi)
{
    CSA *sa = csa_ptr(self);
    i64 i, n;
    int c;

    i = NUM2LONG(oi);
    n = sa->n;

    if (i < 0 || i > n) {    // error
      return Qnil;
    }

    c = sa->BW(sa, i);

    return INT2FIX(c);
}

static VALUE
rb_csa_psi(VALUE self, VALUE oi)
{
    CSA *sa = csa_ptr(self);
    i64 i, n;
    i64 p;

    i = NUM2LONG(oi);
    n = sa->n;

    if (i < 0 || i > n) {    // error
      return Qnil;
    }

    p = sa->psi(sa, i);

    return LONG2FIX(p);
}

static VALUE
rb_csa_lf(VALUE self, VALUE oi)
{
    CSA *sa = csa_ptr(self);
    i64 i, n;
    i64 p;

    i = NUM2LONG(oi);
    n = sa->n;

    if (i < 0 || i > n) {    // error
      return Qnil;
    }

    p = sa->LF(sa, i);

    return LONG2FIX(p);
}

static VALUE
rb_csa_head(VALUE self, VALUE oi)
{
    CSA *sa = csa_ptr(self);
    i64 i, n;
    int c;

    i = NUM2LONG(oi);
    n = sa->n;

    if (i < 0 || i > n) {    // error
      return Qnil;
    }

    c = sa->head(sa, i);

    return INT2FIX(c);
}


static VALUE
rb_csa_lookup(VALUE self, VALUE oi)
{
    CSA *sa = csa_ptr(self);
    i64 i, j, n;

    i = NUM2LONG(oi);
    n = sa->n;

    if (i < 0 || i > n) {    // error
      return Qnil;
    }

    j = sa->lookup(sa, i);

    return LONG2FIX(j);
}

static VALUE
rb_csa_inverse(VALUE self, VALUE oi)
{
    CSA *sa = csa_ptr(self);
    i64 i,j,n;

    i = NUM2LONG(oi);
    n = sa->n;

    if (i < 0 || i > n) {    // error
      return Qnil;
    }

    j = sa->inverse(sa, i);
    return LONG2FIX(j);
}

static VALUE
rb_csa_text(VALUE self, VALUE oi, VALUE oj)
{
    CSA *sa = csa_ptr(self);
    i64 i,j,n;
    uchar *buf;

    i = NUM2LONG(oi);
    j = NUM2LONG(oj);
    n = sa->n;

    if (i < 0 || i > n || j < 0 || j > n) {    // error
      return Qnil;
    }

    buf = (uchar *)alloca(j-i+1+1);
    sa->text(buf, sa, i, j);
    buf[j-i+1] = 0;
    return rb_str_new(buf, j-i+1);
}

static VALUE
rb_csa_substring(VALUE self, VALUE orank, VALUE olen)
{
    CSA *sa = csa_ptr(self);
    i64 rank,len,n;
    uchar *buf;

    rank = NUM2LONG(orank);
    len = NUM2LONG(olen);
    n = sa->n;

    if (rank < 0 || rank > n || len < 0) {    // error
      return Qnil;
    }

    buf = (uchar *)alloca(len+1);
    len = sa->substring(buf, sa, rank, len);
    buf[len] = 0;
    return rb_str_new(buf, len);
}


static VALUE
rb_csa_search(VALUE self, VALUE okey)
{
    CSA *sa = csa_ptr(self);
    i64 n, l, keylen;
    i64 i[2];
    uchar *key;

    key = StringValueCStr(okey);
    keylen = RSTRING(okey)->len;

    if (sa->search(key, keylen, sa, &i[0], &i[1]) < keylen) return Qnil;

    return rb_ary_new3(2, LONG2FIX(i[0]), LONG2FIX(i[1]));
}

static VALUE
rb_csa_search_l(VALUE self, VALUE oc, VALUE range)
{
    CSA *sa = csa_ptr(self);
    i64 ll,rr;
    int c;
    i64 ret;

    c = NUM2INT(oc);
    if (RARRAY(range)->len != 2) {
      return Qnil;
    }
    ll = NUM2LONG(RARRAY(range)->ptr[0]);
    rr = NUM2LONG(RARRAY(range)->ptr[1]);

    ret = sa->searchsub(c, sa, &ll, &rr);
    if (ret == -1) return Qnil;

    return rb_ary_new3(2, LONG2FIX(ll), LONG2FIX(rr));
}

static VALUE
rb_csa_child_l(VALUE self, VALUE range)
{
    CSA *sa = csa_ptr(self);
    i64 l,r,ll,rr;
    int c,i;
    VALUE charset;
    i64 ret;

    if (RARRAY(range)->len != 2) {
      return Qnil;
    }
    ll = NUM2LONG(RARRAY(range)->ptr[0]);
    rr = NUM2LONG(RARRAY(range)->ptr[1]);
    l = r = -1;

    if (!rb_block_given_p()) charset = rb_ary_new();
    for (i=0; i<sa->m; i++) {
      c = sa->AtoC[i];
      l = ll;  r = rr;
      ret = sa->searchsub(c, sa, &l, &r);
      if (ret == 0) {
        if (rb_block_given_p()) {
            rb_yield(rb_ary_new3(2,INT2FIX(c),
                     rb_ary_new3(2,LONG2FIX(l),LONG2FIX(r))));
        } else {
          rb_ary_push(charset,
            rb_ary_new3(2,INT2FIX(c),
            rb_ary_new3(2,LONG2FIX(l),LONG2FIX(r))));
        }
      }
    }
    if (rb_block_given_p()) {
      return Qnil;
    } else {
      return charset;
    }
}

static VALUE
rb_csa_search_r(VALUE self, VALUE olen, VALUE oc, VALUE range)
{
    CSA *sa = csa_ptr(self);
    i64 ll,rr;
    int c;
    i64 len;

    c = NUM2INT(oc);
    if (RARRAY(range)->len != 2) {
      return Qnil;
    }
    ll = NUM2LONG(RARRAY(range)->ptr[0]);
    rr = NUM2LONG(RARRAY(range)->ptr[1]);
    len = NUM2LONG(olen);

    if (csa_search_r(len, c, sa, &ll, &rr) < len) return Qnil;

    return rb_ary_new3(2, LONG2FIX(ll), LONG2FIX(rr));
}

static VALUE
rb_csa_child_r(VALUE self, VALUE olen, VALUE range)
{
    CSA *sa = csa_ptr(self);
    i64 l,r,ll,rr;
    int i,c;
    i64 len;
    VALUE charset;

    if (RARRAY(range)->len != 2) {
      return Qnil;
    }
    ll = NUM2LONG(RARRAY(range)->ptr[0]);
    rr = NUM2LONG(RARRAY(range)->ptr[1]);
    len = NUM2LONG(olen);

    if (!rb_block_given_p()) charset = rb_ary_new();
    for (i=0; i<sa->m; i++) {
      c = sa->AtoC[i];
      l = ll;  r = rr;
      if (csa_search_r(len, c, sa, &l, &r) == len+1) {
        if (rb_block_given_p()) {
          rb_yield(rb_ary_new3(2,INT2FIX(c),
                   rb_ary_new3(2,LONG2FIX(l),LONG2FIX(r))));
        } else {
          rb_ary_push(charset,
                      rb_ary_new3(2,INT2FIX(c),
                      rb_ary_new3(2,LONG2FIX(l),LONG2FIX(r))));
        }
      }
    }
    if (rb_block_given_p()) {
      return Qnil;
    } else {
      return charset;
    }
}

static VALUE
rb_csa_child(VALUE self, VALUE olen, VALUE range, VALUE lr)
{
  if (NUM2INT(lr) == 0) {
    return rb_csa_child_l(self, range);
  } else {
    return rb_csa_child_r(self, olen, range);
  }
}

#if 0
static VALUE
rb_csa_approx_dp(VALUE self, VALUE oc, VALUE od, VALUE okey)
{
    CSA *sa = csa_ptr(self);
    int i;
    int c;
    VALUE *d, *d2;
    int len, keylen;
    uchar *key;
    VALUE ret;

    c = NUM2INT(oc);
    len = RARRAY(od)->len;
    d = RARRAY(od)->ptr;
    key = StringValueCStr(okey);
    keylen = RSTRING(okey)->len;

    d2 = mymalloc(len * sizeof(VALUE));
    d2[0] = INT2FIX(FIX2INT(d[0])+1);

    for (i=1; i<=keylen; i++) {
      d2[i] = INT2FIX(min(FIX2INT(d[i-1]) + ((c != key[keylen-i])?1:0),
                      min(FIX2INT(d[i])+1,
                          FIX2INT(d2[i-1])+1)));
      
    }

    ret = rb_ary_new4(len, d2);
    myfree(d2,len * sizeof(VALUE));
    
    return ret;
}
#endif


void
Init_csa(void)
{
    /* printf("initialize module csa\n"); */
    VALUE rb_cCSA = rb_define_class("CSA", rb_cObject);

    rb_define_alloc_func(rb_cCSA, rb_csa_alloc_func);
    rb_define_method(rb_cCSA, "initialize", rb_csa_initialize, 2);
    rb_define_method(rb_cCSA, "getn", rb_csa_getn, 0);
    rb_define_method(rb_cCSA, "length", rb_csa_getn, 0);
    rb_define_method(rb_cCSA, "size", rb_csa_getn, 0);
    rb_define_method(rb_cCSA, "getm", rb_csa_getm, 0);
    rb_define_method(rb_cCSA, "atoc", rb_csa_atoc, 1);
    rb_define_method(rb_cCSA, "ctoa", rb_csa_ctoa, 1);
    rb_define_method(rb_cCSA, "charset", rb_csa_charset, 0);
    rb_define_method(rb_cCSA, "lookup", rb_csa_lookup, 1);
    rb_define_method(rb_cCSA, "inverse", rb_csa_inverse, 1);
    rb_define_method(rb_cCSA, "text", rb_csa_text, 2);
    rb_define_method(rb_cCSA, "substring", rb_csa_substring, 2);
    rb_define_method(rb_cCSA, "search", rb_csa_search, 1);
    rb_define_method(rb_cCSA, "search_l", rb_csa_search_l, 2);
    rb_define_method(rb_cCSA, "search_r", rb_csa_search_r, 3);
    rb_define_method(rb_cCSA, "child_l", rb_csa_child_l, 1);
    rb_define_method(rb_cCSA, "child_r", rb_csa_child_r, 2);
    rb_define_method(rb_cCSA, "child", rb_csa_child, 3);
//    rb_define_method(rb_cCSA, "approx_dp", rb_csa_approx_dp, 3);
    rb_define_method(rb_cCSA, "bw", rb_csa_bw, 1);
    rb_define_method(rb_cCSA, "psi", rb_csa_psi, 1);
    rb_define_method(rb_cCSA, "lf", rb_csa_lf, 1);
    rb_define_method(rb_cCSA, "head", rb_csa_head, 1);
}
