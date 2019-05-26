
#define ENVINC	64

typedef struct {
 float vco_inc, vco_k;

 float vcf_cutoff, vcf_envmod, vcf_envdecay, vcf_reso, vcf_rescoeff;
 float vcf_e0, vcf_e1;
 float vcf_c0;
 float vcf_d1, vcf_d2;
 float vcf_a, vcf_b, vcf_c;
 int vcf_envpos;

 float vca_attack, vca_decay, vca_a0, vca_a;
 int vca_mode; 
} tb303_t;

#define NOTE (1)
#define A    (1<<1)
#define S    (1<<2)
#define CUT  (1<<3)
#define RES  (1<<4)
#define ENV  (1<<5)
#define DEC  (1<<6)
#define ACC  (1<<7)

typedef struct {
 int note, a, s;
 float cut, res, env, dec, acc;
 int mask;
} tb303_event_t;

void init_tb303(tb303_t *s);
void tb303_lump(float *buf, int len, tb303_t *s);
void tb303_event(tb303_t *s, tb303_event_t *e);
