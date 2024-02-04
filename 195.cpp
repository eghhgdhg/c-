#include<stdio.h>
#include<math.h>
#define PI 3.1416 
#define P0 0.101325     //定义符号常量PI  和大气压力（单位Mpa） 大写以区分         
 float m3=1.38,n=2000,m4=1.81,d=0.095,mk1=0.563,mk2=1.15,rho=0.0415,r=0.0575,l=0.175,m1=0.543;            //定义全局变量 
 int main()
{
	//首先对各个变量以及单位进行说明  质量组（单位为kg） 活塞组质量m3 连杆组质量m4  曲柄销质量mk1 曲柄臂质量mk2  连杆小头质量M1[i] 连杆大头质量m2  旋转质量mr  往复运动质量mj   
	//长度组（单位m） 曲柄臂质心回转半径rho  曲柄半径r  连杆长度l 
	//转速n(转/分钟) 角速度w(rad/s)  曲轴转角alpha[i]  连杆运动角速度wl[i]  
	// 压力组(Mpa)  气缸上部绝对压力p1g[10*i]    大气压力p0  相对压力pg[i]  活塞往复惯性压力pj[i]  活塞上总压力p  侧压力pn[i] 连杆力pl[i]  切向力t 径向力k 
	//径杆比 lambda    
	//主轴颈扭矩(N*m) Mz  连杆轴径扭矩Mq  M为单缸扭矩 
	//连杆轴径负荷 以及方向（Mpa）  pqx[i]   pqy[i] 连杆旋转质量m2的离心力krl   连杆作用于曲柄销上的负荷pq 	alpha[i]q1  pq对应的角度 alpha[i]q   连杆轴承负荷pp  	连杆轴承负荷的对应的角度alphap[i]
	// 主轴颈上的负荷以及方向（Mpa） mk旋转产生的离心力产生的压强krk    主轴颈上的负荷pz    主轴颈上的负荷对应的角度alpha[i]z  pzx[i] pzy[i]  主轴承上的负荷pc   主轴承上的负荷对应的角度alpha[i]c 
	  
	  
	 float alphar[72],alpha[72];
	   float p1g[72]={0.0857,0.08941,0.08336,0.08385,0.08229,0.0818,0.08414,0.08453,
                  0.08346,0.0858,0.08629,0.08751,0.08921,0.09195,0.09239,0.09292,
                  0.09478,0.09582,0.098,0.0986,0.09982,0.10336,0.10736,0.11673,
                  0.13087,0.15526,0.17838,0.21916,0.27409,0.35965,0.49906,0.72989,
                  1.1332,1.87837,3.14011,4.69971,6.70251,6.73661,5.11565,3.51069,
                  2.41675,1.67978,1.2474,0.95628,0.76555,0.64058,0.55072,0.48428,
                  0.44136,0.40536,0.38009,0.36594,0.34233,0.31502,0.20931,0.15575,
                  0.1379,0.11439,0.09486,0.08219,0.09156,0.08953,0.08882,0.09741,
                  0.10317,0.09468,0.08941,0.08424,0.08541,0.09302,0.09673,0.0936};
	  int i;
	 float w,lambda,m2,mk,mr,mj,A;
	  float B[72],beta[72],wl[72],al[72],bet[72],
	        x[72],v[72],a[72], 
	        pl[72],pj[72],p[72],pg[72] ,t[72],k[72] ,pn[72] , 
			Mz1[72],Mz2[72],Mq1[72],M1[72],
			pqx[72] , pqy[72], krl , pq[72] ,	alphaq1[72]  , alphaq[72],pp[72],alphap[72],
			krk,pz[72],alphaz[72],pzx[72],pzy[72],alphaz1[72],pc[72],alphac[72]; 
	                      //定义下变量
	                      
	float sign( float x);//对调用函数进行说明 
	  
	
	  for(i=0;i<=71;i++)
	   {
	   alphar[i]=10*i; 
	   }//对转角的数据进行赋值
	   
	    //一些通用量的计算             
	  w=2*PI*n/60;
	  lambda=r/l;
	  m2=m4-m1;
	  mk=mk1+2*rho/r*mk2;
	  mr=m2+mk;
	  mj=m1+m3;
	  A=PI/4*d*d;
	 
	                    
     for(i=0;i<=71;i++)
	 {printf("曲轴转角=%f  气缸上部绝对压力=%f\n  ",alphar[i],p1g[i]);
	 alpha[i]=alphar[i]*PI/180; 
	 beta[i]=asin(lambda*sin(alpha[i])) ;     //输入公式  第一问   基本上都是顺序结构 
	 bet[i]=180*beta[i]/PI;  
	  B[i]= 1-pow(lambda,2)*pow(sin(alpha[i]),2);                 
     wl[i]=lambda*w*cos(alpha[i])/sqrt(B[i]);                      //所有的三角函数自变量单位是弧度   反三角函数实数即可  
     al[i]=w*w*lambda*(pow(lambda,2)-1)*sin(alpha[i])/pow(B[i],3/2);
    printf(" 摆动角位移(度)beta=%f\n 连杆摆动的角速度(rad/s)wl=%f\n 连杆摆动的角加速度(rad/s*s)al=%f\n",bet[i],wl[i],al[i]);
     
  
	x[i]=r*(1/lambda+1-cos(alpha[i])-1/lambda*cos(beta[i]));    //第二问 
	v[i]=r*w*(sin(alpha[i]+beta[i])/cos(beta[i]));
	a[i]=w*w*r*(cos(alpha[i]+beta[i])/cos(beta[i])+lambda*pow(cos(alpha[i]),2)/pow(cos(beta[i]),3)); 
	 printf(" 活塞位移(m)x=%f\n 活塞速度(m/s)v=%f\n 活塞加速度(m/s*s)a=%f\n",x[i],v[i],a[i]);                                                     
   
    pj[i]=-a[i]*mj/A*pow(10,-6);   //第三问
    pg[i]=p1g[i]-P0;//注意P0是大写 
    p[i]=pg[i]+pj[i];
    pn[i]=p[i]*tan(beta[i]);
    pl[i]=p[i]/cos(beta[i]);
    t[i]=pl[i]*sin(alpha[i]+beta[i]);
    k[i]=pl[i]*cos(alpha[i]+beta[i]);
	printf(" 活塞上总作用力(Mpa)p=%f\n 活塞上侧压力(MPA)pn=%f\n 连杆力(Mpa)pl=%f\n 切向力（Mpa)t=%f\n 径向力（Mpa)k=%f\n",p[i],pn[i],pl[i],t[i],k[i]);        
	
	
	M1[i]=t[i]*r*pow(10,6)*A;                                           //第四问 
	Mz1[i]=0;
	Mz2[i]=M1[i];
	Mq1[i]=Mz1[i]+0.5*M1[i];	
	 printf(" 第一主轴颈扭矩(N*m)Mz1=%f\n 第二主轴颈扭矩(N*m)Mz2=%f\n 连杆轴颈扭矩(N*m)Mq1=%f\n",Mz1[i],Mz2[i],Mq1[i]);
												          

														   //第五问 
	
	krl=m2*r*w*w/A*pow(10,-6);													   
	pqx[i]=krl-k[i];
	pqy[i]=t[i];													   																									   
	pq[i]=sqrt(pqx[i]*pqx[i]+pqy[i]*pqy[i]);
	alphaq1[i]=atan(fabs(pqy[i]/pqx[i]));
	alphaq[i]=(1-sign(pqx[i]))*PI/2+sign(pqx[i])*sign(pqy[i])*alphaq1[i];		
	pp[i]=pq[i];
	alphap[i]=alphaq[i]+PI+alpha[i]+beta[i];									           
	
	krk=mk*r*w*w/A*pow(10,-6); //由于只是单拐，所以pz1=pz2,这里只计算一个，写为pz 
	pzx[i]=-0.5*(krl+krk-k[i]);
	pzy[i]=-0.5*t[i];
	pz[i]=sqrt(pzx[i]*pzx[i]+pzy[i]*pzy[i]);
	alphaz1[i]=atan(fabs(pzy[i]/pzx[i]));
	alphaz[i]=(1-sign(pzx[i]))*PI/2+sign(pzx[i])*sign(pzy[i])*alphaz1[i];
	pc[i]=pz[i];
	alphac[i]=alphaz[i]+PI+alpha[i];
	printf(" 主轴颈上的负荷(Mpa)pz=%f\n 主轴颈上的负荷对应的角度(rad)alphaz=%f\n 主轴承上的负荷(Mpa)pc=%f\n  主轴承上的负荷对应的角度（rad)alphac=%f\n ",pz[i],alphaz[i],pc[i],alphac[i]);
    printf(" 连杆轴颈上的负荷(Mpa)pq=%f\n 连杆轴颈上的负荷对应的角度(rad)alphaq=%f\n 连杆轴承上的负荷(Mpa)pp=%f\n  连杆轴承上的负荷对应的角度（rad)alphap=%f\n ",pq[i],alphaq[i],pp[i],alphap[i]);
    printf("\n\n");
    
	 
}
//printf("%f  %f\n",pj[220],pg[220]);
//最后需要一个连杆轴径轴承，主轴颈轴承的负荷矢量图
for(i=0;i<=71;i++)
	 {alphaz[i]=180*alphaz[i]/PI;
	 alphac[i]=180*alphac[i]/PI;
	 alphaq[i]=180*alphaq[i]/PI;
	 alphap[i]=180*alphap[i]/PI;
	 p[i]=p[i]*pow(10,6)*A;
	 /*pz[i]=pz[i]*pow(10,6)*A;
	  pc[i]=pc[i]*pow(10,6)*A;
	   pq[i]=pq[i]*pow(10,6)*A;
	    pp[i]=pp[i]*pow(10,6)*A;*/
	    p[i]=p[i]*pow(10,6)*A; 
	 printf("%12f %12f %12f %12f %12f %12f %12f %12f %12f %12f %12f\n",alphar[i],p1g[i],pz[i],alphaz[i],pc[i],alphac[i],pq[i],alphaq[i],pp[i],alphap[i],p[i]);
}   					           
	return 0; 
}
 
 
 float sign( float x)
 { float z;
 z=x/fabs(x);
 return z;
 }

