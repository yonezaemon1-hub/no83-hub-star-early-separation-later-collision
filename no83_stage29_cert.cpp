#include <cstdio>
#include <array>
#include <vector>
extern "C" {
typedef long mpfr_prec_t; typedef long mpfr_sign_t; typedef long mpfr_exp_t;
typedef struct {mpfr_prec_t _mpfr_prec; mpfr_sign_t _mpfr_sign; mpfr_exp_t _mpfr_exp; unsigned long *_mpfr_d;} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1];
int mpfr_div(mpfr_t,const mpfr_t,const mpfr_t,int);void mpfr_init2(mpfr_t,long);void mpfr_clear(mpfr_t);int mpfr_set(mpfr_t,const mpfr_t,int);int mpfr_set_ui(mpfr_t,unsigned long,int);int mpfr_div_ui(mpfr_t,const mpfr_t,unsigned long,int);int mpfr_sqrt(mpfr_t,const mpfr_t,int);int mpfr_asin(mpfr_t,const mpfr_t,int);int mpfr_const_pi(mpfr_t,int);int mpfr_sub(mpfr_t,const mpfr_t,const mpfr_t,int);int mpfr_add(mpfr_t,const mpfr_t,const mpfr_t,int);int mpfr_mul(mpfr_t,const mpfr_t,const mpfr_t,int);int mpfr_sin(mpfr_t,const mpfr_t,int);int mpfr_neg(mpfr_t,const mpfr_t,int);double mpfr_get_d(const mpfr_t,int);
}
constexpr int DN=3,UP=2;
struct I{mpfr_t l,u;I(){mpfr_init2(l,160);mpfr_init2(u,160);}~I(){mpfr_clear(l);mpfr_clear(u);}I(const I&)=delete;I&operator=(const I&)=delete;};
void cp(I&a,const I&b){mpfr_set(a.l,b.l,DN);mpfr_set(a.u,b.u,UP);}void ui(I&a,unsigned long n){mpfr_set_ui(a.l,n,DN);mpfr_set_ui(a.u,n,UP);}void neg(I&a,const I&b){mpfr_neg(a.l,b.u,DN);mpfr_neg(a.u,b.l,UP);}void add(I&a,const I&b,const I&c){mpfr_add(a.l,b.l,c.l,DN);mpfr_add(a.u,b.u,c.u,UP);}void sub(I&a,const I&b,const I&c){mpfr_sub(a.l,b.l,c.u,DN);mpfr_sub(a.u,b.u,c.l,UP);}void divui(I&a,const I&b,unsigned long n){mpfr_div_ui(a.l,b.l,n,DN);mpfr_div_ui(a.u,b.u,n,UP);}
void sn(I&a,const I&b){I width;mpfr_sub(width.u,b.u,b.l,UP);mpfr_sin(a.l,b.l,DN);mpfr_sin(a.u,b.l,UP);mpfr_sub(a.l,a.l,width.u,DN);mpfr_add(a.u,a.u,width.u,UP);}
void init(std::array<I,5>&x){I one,three,r,alpha,pi,tmp;ui(one,1);ui(three,3);mpfr_sqrt(r.l,three.l,DN);mpfr_sqrt(r.u,three.u,UP);mpfr_div(tmp.l,one.l,r.u,DN);mpfr_div(tmp.u,one.u,r.l,UP);mpfr_asin(alpha.l,tmp.l,DN);mpfr_asin(alpha.u,tmp.u,UP);mpfr_const_pi(pi.l,DN);mpfr_const_pi(pi.u,UP);ui(x[0],0);cp(x[1],alpha);sub(x[2],pi,alpha);neg(x[3],alpha);neg(x[4],alpha);}
using V=std::array<I,5>;using E=std::vector<std::array<int,2>>;
void deriv(V&f,const V&x,const E&edges){for(auto&v:f)ui(v,0);I diff,s,minus,tmp;for(auto e:edges){sub(diff,x[e[1]],x[e[0]]);sn(s,diff);add(tmp,f[e[0]],s);cp(f[e[0]],tmp);neg(minus,s);add(tmp,f[e[1]],minus);cp(f[e[1]],tmp);}for(auto&v:f){divui(tmp,v,5);cp(v,tmp);}}
void step(V&x,const E&edges){V f,mid,g;I tmp,half,delta;deriv(f,x,edges);for(int i=0;i<5;i++){divui(half,f[i],20000);add(tmp,x[i],half);cp(mid[i],tmp);}deriv(g,mid,edges);for(int i=0;i<5;i++){divui(delta,g[i],10000);add(tmp,x[i],delta);cp(x[i],tmp);}}
int main(){V x,y;init(x);init(y);E hub={{0,1},{0,2},{0,3},{0,4}};E e0=hub,e1=hub;for(auto e:E{{1,3},{1,4},{2,3},{2,4}})e0.push_back(e);e1.push_back({1,2});I d;for(int k=1;k<=31200;k++){step(x,e0);step(y,e1);if(k==31000||k==31200){sub(d,y[0],x[0]);printf("step %d D_midpoint interval [%.17g, %.17g]\n",k,mpfr_get_d(d.l,DN),mpfr_get_d(d.u,UP));fflush(stdout);}}}
