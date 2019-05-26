#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tb303.h"

void init_tb303(tb303_t *s) {
 s->vco_inc = 0.0f;
 s->vco_k = 0.0f;

 s->vcf_cutoff =
  s->vcf_envmod = 
  s->vcf_reso = 
  s->vcf_envdecay = 
  s->vcf_a =
  s->vcf_b = 
  s->vcf_d1 =
  s->vcf_d2 = 
  s->vcf_c0 =
  s->vcf_e0 =
  s->vcf_e1 =
  s->vcf_a = 0.0f;

 s->vca_mode = 2;
 s->vca_attack = 1.0f - 0.94406088f;
 s->vca_decay = 0.99897516;
 s->vca_a0 = 0.5f;
}

void tb303_event(tb303_t *s, tb303_event_t *e) {
 int dirty = 0;
 float d;

 if( (e->mask & CUT) == CUT) {
  s->vcf_cutoff = e->cut;
  dirty = 1;
 }
 if( (e->mask & RES) == RES) {
  s->vcf_reso = e->res;
  s->vcf_rescoeff = exp(-1.20 + 3.455*s->vcf_reso);
  dirty = 1; 
 }
 if( (e->mask & ENV) == ENV) {
  s->vcf_envmod = e->env;
  dirty = 1;
 }
 if( (e->mask & DEC) == DEC) {
  d = e->dec;
  d = 0.2 + (2.3*d);
  d*=44100;
  s->vcf_envdecay = pow(0.1, 1.0/d * ENVINC);
 }
 if(dirty == 1) {
  s->vcf_e1 = exp(6.109 + 1.5876* s->vcf_envmod + 
                  2.1553*s->vcf_cutoff - 1.2*(1.0-s->vcf_reso));
  s->vcf_e0 = exp(5.613 - 0.8*s->vcf_envmod + 
                  2.1553*s->vcf_cutoff - 0.7696*(1.0-s->vcf_reso));
  s->vcf_e0*=M_PI/44100.0f;
  s->vcf_e1*=M_PI/44100.0f;
  s->vcf_e1 -= s->vcf_e0;
  s->vcf_envpos = ENVINC;
 }
 if( (e->mask & NOTE) == NOTE) {
  s->vco_inc = (440.0/44100.0)*pow(2, (e->note-57)*(1.0/12.0));
  s->vca_mode = 0;
  s->vcf_c0 = s->vcf_e1;
  s->vcf_envpos = ENVINC; 
 } else
  s->vca_mode = 1;
}

void tb303_lump(float *buf, int len, tb303_t *s) {
 int i;
 float w,k;
 for(i=0;i<len;i++) {

  if(s->vcf_envpos >= ENVINC) {
   w = s->vcf_e0 + s->vcf_c0;
   k = exp(-w/s->vcf_rescoeff);
   s->vcf_c0 *= s->vcf_envdecay;
   s->vcf_a = 2.0*cos(2.0*w) * k;
   s->vcf_b = -k*k;
   s->vcf_c = 1.0 - s->vcf_a - s->vcf_b;
   s->vcf_envpos = 0;
  }
		
  buf[i]=s->vcf_a*s->vcf_d1 + 
   s->vcf_b*s->vcf_d2 + s->vcf_c*s->vco_k*s->vca_a;
  s->vcf_d2=s->vcf_d1; 
  s->vcf_envpos++;
  s->vcf_d1=buf[i];

  s->vco_k += s->vco_inc;
  if(s->vco_k > 0.5) s->vco_k -= 1.0;

  if(i==len/2) s->vca_mode = 2;
  if(!s->vca_mode) 
   s->vca_a+=(s->vca_a0-s->vca_a)*s->vca_attack;
  else 
   if(s->vca_mode == 1) {
    s->vca_a *= s->vca_decay;
    if(s->vca_a < (1/65536.0)) { s->vca_a = 0; s->vca_mode = 2; }
   }

 }
}


