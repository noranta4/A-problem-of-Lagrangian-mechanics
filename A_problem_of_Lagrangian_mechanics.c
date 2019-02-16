/*ANTONIO NORELLI PROF. RICCI TERSENGHI*/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define A 4./(3*M_PI)

/*STRUCT***************************************************/
typedef struct{
	double x, th, vx, vth, dt;
	}crd;
	
typedef struct{
	double md, mp, r, k, g;
	}cost;
	
typedef struct{
	double pot, cin, tot;
	}ene;
		
/*FUNZIONI************************************************/
crd rkstep(crd, cost);
inline double ax(crd, cost);
inline double ath(crd, cost);
inline ene e(crd, cost);
inline lambda(cost);
int peq(crd*);
cost setcost(char*);
crd setci(char*);

int main(){
	srand48(time(0));
	/*SETTAGGIO NOME DEL FILE, COSTANTI, CONDIZIONI INIZIALI*/
	char nfile[14]="anCOST#CI#.txt", f;
	cost j = setcost(&f);
	nfile[6]=f;
	crd p0 = setci(&f);
	nfile[9]=f;
	/*DICHIARAZIONE VARIABILI, TEMPO TOTALE DI INTEGRAZIONE DA INPUT*/
	crd tmp=p0, pre, crdkmax;
	ene set;
	double tist, tmax, k2dt=0, k1dt=0, k0dt=0, kmax=0, eini=e(p0,j).tot;
	int maxfile=1, i=0;
	printf("\nDigitare il tempo totale di integrazione\n");
	scanf("%lf", &tmax);
	/*APERTURA FILE, CONDIZIONE PER CONTENERNE LE DIMENSIONI, INIZIALIZZAZIONE*/
	FILE *fp;
	if(tmax>10000*p0.dt){
		maxfile=(tmax/p0.dt)/10000;
	}
	if((fp=fopen(nfile,"w")) == NULL){
		printf("errore nell'apertura del file");
	}
	fprintf(fp,"#File di dati\n#md=%.4lf mp=%.4lf r=%.4lf k=%.4lf g=%.4lf dt=%lf tmax=%lf \n#1:t 2:x 3:theta 4:X 5:Y 6:vx 7:vth 8:U 9:K 10:E 11:(E(t)-E(0))/E(0)\n", j.md, j.mp, j.r, j.k, j.g, p0.dt, tmax);
	
	for(i=0;i<21;i++){
		if(i!=0){
			tmax=5;
			tmp.dt=0.000005;
		}
		k2dt=k1dt=k0dt=kmax=0;
		tist=0;
		while(tist<=tmax){
			set=e(tmp, j);
			/*CONDIZIONE PER CONTENERE LE DIMENSIONI DEL FILE*/
			if(((int)(tist/p0.dt)%maxfile==0)&&(i==0)){
				fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %.15lf %.15lf\n", tist, tmp.x, tmp.th, tmp.x*cos(tmp.th), tmp.x*sin(tmp.th), tmp.vx, tmp.vth, set.pot, set.cin, set.tot, (set.tot-eini)/eini);
			}
			/*RICERCA MASSIMO ENERGIA CINETICA*/
			k2dt=k1dt;
			k1dt=k0dt;
			k0dt=set.cin;
			if((((k1dt>k2dt)&&(k1dt>k0dt))||(tist==tmax))&&(k1dt>kmax)){
				kmax=k1dt;
				crdkmax=pre;
//printf("\nk=%.20lf %lf", kmax, tist);
				if(i!=0){
					tist=tmax;
				}
			}
			pre=tmp;
			/*PASSO RUNGE-KUTTA2*/
			tmp=rkstep(tmp, j);
			tist+=tmp.dt;
		}
		tmp=crdkmax;
		tmp.vx=0;
		tmp.vth=0;
//printf("\nKM=%.20lf", kmax);
	}
	printf("\nPunto di equilibrio stabile in (theta in rad): \nx=    %.30lf \ntheta=%.30lf\n", tmp.x, tmp.th);
	fclose(fp);
}

/*DOPPIO PASSO RUNGE-KUTTA 2*************************************/
crd rkstep(crd h, cost j){
	crd plus;
	plus.x=h.x+0.5*h.vx*h.dt;
	plus.th=h.th+0.5*h.vth*h.dt;
	plus.vx=h.vx+0.5*ax(h, j)*h.dt;
	plus.vth=h.vth+0.5*ath(h, j)*h.dt;
	plus.dt=h.dt;
	
	h.x=h.x+plus.vx*h.dt;
	h.th=h.th+plus.vth*h.dt;
	h.vx=h.vx+ax(plus, j)*h.dt;
	h.vth=h.vth+ath(plus, j)*h.dt;
	return h;
}

/*ACCELERAZIONE X E THETA****************************************/
inline double ax(crd h, cost j){
	return h.x*h.vth*h.vth-h.x*(j.k/j.mp)-j.g*sin(h.th);
}
inline double ath(crd h, cost j){
	return (-2*j.mp*h.x*h.vx*h.vth-j.mp*j.g*h.x*cos(h.th)-j.md*j.g*j.r*A*sin(h.th))/(j.mp*h.x*h.x+0.5*j.md*j.r*j.r);
}

/*ENERGIE U, K, TOT*********************************************/
inline ene e(crd h, cost j){
	ene a;
	a.pot=0.5*j.k*h.x*h.x+j.mp*j.g*h.x*sin(h.th)-A*j.r*j.md*j.g*cos(h.th);
	a.cin=0.5*j.mp*h.vx*h.vx+0.5*j.mp*h.x*h.x*h.vth*h.vth+0.25*j.md*j.r*j.r*h.vth*h.vth;
	a.tot=a.pot+a.cin;
	return a;
}

/*GRANDEZZA LAMBDA***********************************************/
inline lambda(cost a){
	return (A*a.md*a.r*a.k)/(a.mp*a.mp*a.g);
}

/*SCELTA COSTANTI************************************************/
cost setcost(char *b){
	do{
		printf("\nSi scelgano le costanti del sistema, digitare:\na   per (lambda>1) md=10, mp=1, k=1, r=8, g=3*pi\nb   per (lambda<1) md=1, mp=1, k=1, r=8, g=3*pi\nc   per sceglierle a piacere\nd   per sceglierle a caso\n\n");
		scanf("%c%c", b);
	}
	while((*b<97)||(*b>100));
	if(*b==97){
		cost j={10, 1, 8, 1, 3*M_PI};
		return j;
	}
	else if(*b==98){
		cost j={1, 1, 8, 1, 3*M_PI};
		return j;
	}
	else if (*b==99){
		cost j;
		printf("\nmd= ");
		scanf("%lf", &j.md);
		printf("\nmp= ");
		scanf("%lf", &j.mp);
		printf("\nk= ");
		scanf("%lf", &j.k);
		printf("\nr= ");
		scanf("%lf", &j.r);
		printf("\ng= ");
		scanf("%lf%c", &j.g);
		return j;
	}
	else{
		cost j={100*drand48(), 100*drand48(), 100*drand48(), 100*drand48(), 100*drand48()};
		printf("Le costanti scelte sono:\nmd=%.4lf mp=%.4lf r=%.4lf k=%.4lf g=%.4lf\n", j.md, j.mp, j.r, j.k, j.g);
		if(lambda(j)>1){
			printf("(lambda>1)\n");
		}
		else{
			printf("(lambda<=1)\n");
		}
		return j;
	}
		
	
}

/*SCELTA CONDIZIONI INIZIALI***********************************************************/	
crd setci(char *c){
	if(*c==97){
		do{
			printf("\nSi scelgano le condizioni iniziali, digitare:\na   per sceglierle a piacere\nb   per scegliere una config. con oscillazioni in fase\nc   per scegliere una config. con oscillazioni in controfase\nd   per sceglierle a caso\n\n");
			scanf("%c%c", c);
		}
		while((*c<97)||(*c>100));
		if(*c==97){
			crd ci;
			printf("\nx0= ");
			scanf("%lf", &ci.x);
			printf("\ntheta0= ");
			scanf("%lf", &ci.th);
			printf("\nvx0= ");
			scanf("%lf", &ci.vx);
			printf("\nvtheta0= ");
			scanf("%lf", &ci.vth);
			printf("\n\nDigitare dt (si consiglia 0.000005)\n");
			scanf("%lf", &ci.dt);
			return ci;
		}
		else if(*c==98){
			crd ci={10.1, 0.5, 0, 0, 0};
			printf("\n\nDigitare dt (si consiglia 0.000005)\n");
			scanf("%lf", &ci.dt);
			return ci;
		}
		else if(*c==99){
			crd ci={-9.3, 0.5, 0, 0, 0};
			printf("\n\nDigitare dt (si consiglia 0.000005)\n");
			scanf("%lf", &ci.dt);
			return ci;
		}
		else{
			crd ci={100*drand48(), 2*M_PI*drand48(), 100*drand48(), 100*drand48(), 0};
			printf("Le cond. in. scelte sono:\nx=%.4lf theta=%.4lf vx=%.4lf vtheta=%.4lf\n", ci.x, ci.th, ci.vx, ci.vth);
			printf("\nDigitare dt (si consiglia 0.000005)\n");
			scanf("%lf", &ci.dt);
			return ci;
		}
	}
	else{
		do{
			printf("\nSi scelgano le condizioni iniziali, digitare:\na   per sceglierle a piacere\nd   per sceglierle a caso\n\n");
			scanf("%c%c", c);
		}
		while((*c!=97)&&(*c!=100));
		if(*c==97){
			crd ci;
			printf("\nx0= ");
			scanf("%lf", &ci.x);
			printf("\ntheta0= ");
			scanf("%lf", &ci.th);
			printf("\nvx0= ");
			scanf("%lf", &ci.vx);
			printf("\nvtheta0= ");
			scanf("%lf", &ci.vth);
			printf("\n\nDigitare dt (si consiglia 0.000005)\n");
			scanf("%lf", &ci.dt);
			return ci;
		}
		else{
			crd ci={1000*drand48(), 2*M_PI*drand48(), 1000*drand48(), 1000*drand48(), 0};
			printf("Le cond. in. scelte sono:\nx=%.4lf theta=%.4lf vx=%.4lf vtheta=%.4lf\n", ci.x, ci.th, ci.vx, ci.vth);
			printf("\nDigitare dt (si consiglia 0.000005)\n");
			scanf("%lf", &ci.dt);
			return ci;
		}
	}
}
