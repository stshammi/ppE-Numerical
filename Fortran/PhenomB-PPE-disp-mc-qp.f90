program bd_lisa_noangle

use nrtype_qp
implicit none


integer, parameter :: imax=7, nbin=1000,nmass=100,ndet=2,nb1=3,nb2=8, &
        nfmax=1 !1:f3, 2:132Hz, 3:ISCO
integer :: k,l,m,n
real(dp) :: m1, m2, mtotal, mc, eta, tsolar=4.925491_dp*10.0_dp**(-6), &
        dl, beta, tc, phic, tyear=3.1536_dp*10.0_dp**7, sigma, &
        fmax, fmin, tobs=1.0_dp, s=0.3_dp, &
        fend, fcut, fisco, fin,c=2.9979_dp*10.0_dp**10, &
        z, betag, d, snthr, integ, amp, kpc, dsub, sn,chi,f,f3,noiseligo, &
        lambdag,const,SNR_des,b,nPNmin,fend_insp,factor,epsilon
real(dp), dimension(nbin) :: w, fbin
real(dp), dimension(imax,imax) :: fm, fmi!,check,fmnorm,fminorm
complex(dpc) :: ad1,ad2,add1,add2
!real(dp), dimension(imax) :: diag

!complex(dpc) :: ad1

open (unit=1,file="Data/PhenomB-PPE-full-20Hz-nbin1000-eps20-f3-disp-mc-qp.dat")

!sigma90=1.64485

SNR_des=24._dp

!fend=10.0_dp**(4.0_dp)
fcut=20.0_dp**(1.0_dp)
fend_insp=132.0_dp
z=0.088_dp

chi=0.0_dp

!do k=0,nmass

!snthr=5.0_dp

!m1=10.0_dp**(6.0_dp+3.0_dp/nmass*k)*tsolar
!m1=10.0_dp**(4.0_dp)*tsolar

!m1=10.0_dp**3_dp*tsolar
!m2=10.0_dp**4
!m2=m1

!mtotal=m1+m2
!eta=m1*m2/mtotal**2

!mtotal=10.0_dp**(1.0_dp+2.0_dp/nmass*k)*tsolar

m1=35.7_dp*(1.0_dp+z)
m2=29.1_dp*(1.0_dp+z)
mtotal=(m1+m2)*tsolar
eta=m1*m2/(mtotal/tsolar)**2
mc=eta**(3.0_dp/5.0_dp)*mtotal

betag=0.0_dp

! S/N=10
!dl=1.3037_dp*10.0_dp**10*tyear
dl=0.41_dp*10.0_dp**9*3.261636_dp*tyear
!d=1.99299_dp*10.0_dp**9*3.261636_dp*tyear
d=0.377398_dp*10.0_dp**9*3.261636_dp*tyear

!dl=1.0_dp*10.0_dp**4*3.261636_dp*tyear
!d=dl

!d=2.77054_dp*10.0_dp**9*3.261636_dp*tyear
!d=2.0437_dp*10.0_dp**9*3.261636_dp*tyear

!beta=0.0_dp
tc=0.0_dp
phic=0.0_dp
!omegab=0.0_dp
!sigma=0.0_dp
!betag=0.0_dp

factor=SNR_des/51.830158373945181_dp

const=(mtotal*eta**(3.0_dp/5.0_dp))**(5.0_dp/6.0_dp) &
       /(pi_d**(2.0_dp/3.0_dp))/dl*sqrt(5.0_dp/24.0_dp)*factor

!check=0.0_dp

fisco=(6.0_dp**(1.5_dp)*pi_d*mtotal)**(-1)
!fin=4.149_dp*10.0_dp**(-5.0_dp)*(Mc/tsolar/10.0_dp**6)**(-5.0_dp/8.0_dp) &
!     *(tobs)**(-3.0_dp/8.0_dp)

!fmax=min(f3, fend)
!fmin=max(fin, fcut)
!fmax=2.199_dp*10.0_dp**(-4)


do k=0,nb1-1

!nPNmin=-4.0_dp

!b=(nPNmin+k/2.0_dp)*2.0_dp/3.0_dp-5.0_dp/3.0_dp

b=(k/2.0_dp)-1.0_dp

call wave(f,mc,eta,chi,const,tc,phic,betag,ad1,add1,f3)

!fmax=f3
!fmax=132.0_dp
!fmax=fisco

if (nfmax==1) then

fmax=f3
epsilon=10.0_dp**(-20)

else if (nfmax==2) then

fmax=fend_insp
epsilon=10.0_dp**(-9)

else if (nfmax==3) then

fmax=fisco
epsilon=10.0_dp**(-9)

end if


fmin=fcut



print *, "alpha=", k/2.0_dp


call fish(mtotal,eta,chi,const,tc,phic,betag,fm,fmi)

!lambdag=sqrt(pi_d*d*c**2/((1.0_dp+z)*(sqrt(fmi(7,7))*(pi_d*mc)**b)))


write (1,'(2E22.14)') k/2.0_dp,sqrt(fmi(7,7))
!                    
!print *, "nPN=", nPNmin+k/2.0_dp

end do



k=3

b=(k/2.0_dp)-1.0_dp

call wave(f,mc,eta,chi,const,tc,phic,betag,ad1,add1,f3)

!fmax=f3
!fmax=132.0_dp
!fmax=fisco

if (nfmax==1) then

fmax=f3
epsilon=10.0_dp**(-20)

else if (nfmax==2) then

fmax=fend_insp
epsilon=10.0_dp**(-9)

else if (nfmax==3) then

fmax=fisco
epsilon=10.0_dp**(-9)

end if


fmin=fcut



print *, "alpha=", k/2.0_dp


call fish(mtotal,eta,chi,const,tc,phic,betag,fm,fmi)

!lambdag=sqrt(pi_d*d*c**2/((1.0_dp+z)*(sqrt(fmi(7,7))*(pi_d*mc)**b)))


write (1,'(2E22.14)') k/2.0_dp,sqrt(fmi(7,7))
!






do k=nb1+2,nb2

nPNmin=-4.0_dp

!b=(nPNmin+k/2.0_dp)*2.0_dp/3.0_dp-5.0_dp/3.0_dp
b=(k/2.0_dp)-1.0_dp

call wave(f,mc,eta,chi,const,tc,phic,betag,ad1,add1,f3)

!fmax=f3
!fmax=132.0_dp
!fmax=fisco

if (nfmax==1) then

fmax=f3
epsilon=10.0_dp**(-20)

else if (nfmax==2) then

fmax=fend_insp
epsilon=10.0_dp**(-9)

else if (nfmax==3) then

fmax=fisco
epsilon=10.0_dp**(-6)

end if


fmin=fcut


call fish(mtotal,eta,chi,const,tc,phic,betag,fm,fmi)

!lambdag=sqrt(pi_d*d*c**2/((1.0_dp+z)*(sqrt(fmi(7,7))*(pi_d*mc)**b)))


write (1,'(2E22.14)') k/2.0_dp,sqrt(fmi(7,7))
!
print *, "nPN=", nPNmin+k/2.0_dp

end do


print *, "S/N =",sqrt(fm(4,4))
print *, "fmax=",fmax
print *, "fmin=",fmin
!print *, "dl(kpc)=",dl

!print *, "omega_BD(uncor)=", sqrt(fm(7,7))

!print *, fm(1,1)
!print *, fm(2,2)
!print *, fm(1,2)
!print *, fm(2,1)
!print *, fm(3,3)
!print *, fm(4,4)
!print *, fm(5,5)
!print *, fm(6,6)
print *, "delta mtotal/mtotal=", sqrt(fmi(1,1))
print *, "delta eta/eta=", sqrt(fmi(2,2))
print *, "delta chi=", sqrt(fmi(3,3))
print *, "delta dl/dl=", sqrt(fmi(4,4))
print *, "delta tc=", sqrt(fmi(5,5))
print *, "delta phic=", sqrt(fmi(6,6))
print *, "beta_g=", sqrt(fmi(7,7))
!print *, "lambda_g=", lambdag



close(1)

contains



subroutine fish(mtotal,eta,chi,const,tc,phic,betag,fm,fmi)

use nrtype_qp
implicit none

real(dp), intent(in) :: mtotal,eta,chi,const,tc,phic,betag
real(dp), dimension(imax,imax), intent(inout) :: fm,fmi

integer :: i,j,n
real(dp) :: noise, amp, x, f, v, k4, a4, b4, c4, k5, a5, b5, c5, &
         snsa, s_gal, s_exgal, kappa=4.5_dp, dndf, shot, rad, accel, nn, &
         f0, clean=0.01_dp, noisedec, &
         aa0,aa1,aa2,aa3,aa4,aa5,aa6
complex(dpc) :: wi=(0.0_dp, 1.0_dp), adt1, adt2, ad1, ad2, adt_1, adt_2, &
                add1, add2, addt1, addt2, addt_1, addt_2
complex(dpc), dimension(imax) ::dad1, dad2
real(dp), dimension(0,0) :: b
 

! frequency integration

fm=0.0_dp
fmi=0.0_dp
sn=0.0_dp

call gauleg(fmin,fmax,fbin,w,nbin)

do n=1, nbin

f=fbin(n)


call wave(f,mc,eta,chi,const,tc,phic,betag,ad1,add1,f3)

call wave(f,mc*(1._dp+epsilon),eta,chi,const,tc,phic,betag,adt1,addt1,f3)
call wave(f,mc*(1._dp-epsilon),eta,chi,const,tc,phic,betag,adt_1,addt_1,f3)
dad1(1)=(adt1-adt_1)/(2.0_dp*epsilon)!*mtotal
!dad2(1)=(adt2-adt_2)/(2.0_dp*epsilon)

call wave(f,mc,eta*(1._dp+epsilon),chi,const,tc,phic,betag,adt1,addt1,f3)
call wave(f,mc,eta*(1._dp-epsilon),chi,const,tc,phic,betag,adt_1,addt_1,f3)
dad1(2)=(adt1-adt_1)/(2.0_dp*epsilon)!*eta
!dad2(2)=(adt2-adt_2)/(2.0_dp*epsilon)

call wave(f,mc,eta,chi+epsilon,const,tc,phic,betag,adt1,addt1,f3)
call wave(f,mc,eta,chi-epsilon,const,tc,phic,betag,adt_1,addt_1,f3)
dad1(3)=(adt1-adt_1)/(2.0_dp*epsilon)
!dad2(3)=(adt2-adt_2)/(2.0_dp*epsilon)

call wave(f,mc,eta,chi,const*(1._dp+epsilon),tc,phic,betag,adt1,addt1,f3)
call wave(f,mc,eta,chi,const*(1._dp-epsilon),tc,phic,betag,adt_1,addt_1,f3)
dad1(4)=(adt1-adt_1)/(2.0_dp*epsilon)
!dad2(4)=(adt2-adt_2)/(2.0_dp*epsilon)*dl
!write (1,'(5E22.14)') real(dad1(4)),real(adt1),real(adt_1),real(adt1-adt_1),real(ad1)


call wave(f,mc,eta,chi,const,tc+epsilon,phic,betag,adt1,addt1,f3)
call wave(f,mc,eta,chi,const,tc-epsilon,phic,betag,adt_1,addt_1,f3)
dad1(5)=(adt1-adt_1)/(2.0_dp*epsilon)
!dad2(5)=(adt2-adt_2)/(2.0_dp*epsilon)

call wave(f,mc,eta,chi,const,tc,phic+epsilon,betag,adt1,addt1,f3)
call wave(f,mc,eta,chi,const,tc,phic-epsilon,betag,adt_1,addt_1,f3)
dad1(6)=(adt1-adt_1)/(2.0_dp*epsilon)
!dad2(6)=(adt2-adt_2)/(2.0_dp*epsilon)

call wave(f,mc,eta,chi,const,tc,phic,betag+epsilon,adt1,addt1,f3)
call wave(f,mc,eta,chi,const,tc,phic,betag-epsilon,adt_1,addt_1,f3)
dad1(7)=(adt1-adt_1)/(2.0_dp*epsilon)

!call wave2(f,mc,eta,chi,const,tc,phic,betag,adt1,addt1,f3)

! noise curve

!f0=7.36_dp

!shot=2.3_dp*10.0_dp**(-24)*sqrt(1.0_dp+(f/f0)**2)
!rad=6.0_dp*10.0_dp**(-26)*f**(-2)/sqrt(1.0_dp+(f/f0)**2)
!accel=2.0_dp*10.0_dp**(-26)*f**(-2)
!snsa=shot**2+rad**2+accel**2
!nn=clean*1.3_dp*10.0_dp**(-48)*f**(-7.0_dp/3.0_dp)
!gal=2.1_dp*10.0_dp**(-45)*f**(-7.0_dp/3.0_dp)*exp(-(f/0.05_dp)**2)**2
!exgal=4.2_dp*10.0_dp**(-47)*f**(-7.0_dp/3.0_dp)*exp(-(f/0.05_dp)**2)**2
!dndf=2.0_dp*10.0_dp**(-3)*f**(-11.0_dp/3.0_dp)

!*******
!LISA
!*******

!snsa=9.18_dp*10.0_dp**(-52)*f**(-4)+1.59_dp*10.0_dp**(-41) &
!     +9.18_dp*10.0_dp**(-38)*f**2
!s_gal=2.1_dp*10.0_dp**(-45)*f**(-7.0_dp/3.0_dp)
!s_exgal=4.2_dp*10.0_dp**(-47)*f**(-7.0_dp/3.0_dp)
!dndf=2.0_dp*10.0_dp**(-3)*f**(-11.0_dp/3.0_dp)

!noise=min(snsa/exp(-kappa/tyear*dndf),snsa+s_gal)+s_exgal

!print *, dad1(1)


!*******
!DECIGO
!*******

!noisedec=(min(snsa/exp(-kappa*dndf/tyear),snsa+gal)+exgal+nn)*10.0_dp**(16)


!if (f < 5.0_dp*10.0_dp**(-1)) then
!        noise = (10.0_dp**(-15.51_dp)*f**(-1.7_dp))**2
        
!  else  
!        noise = (10.0_dp**(-14.85_dp)*f**(0.5_dp))**2
  
!end if

!*******
!LIGO
!*******

!Ajith et al.

!f0=215.0_dp

!noiseligo=10.0_dp**(-49._dp)*((f/f0)**(-4.14_dp)-5.0_dp*(f/f0)**(-2._dp) &
!            +111.0_dp*(1.0_dp-(f/f0)**2._dp+(f/f0)**4._dp/2._dp) &
!           /(1.0_dp+(f/f0)**2._dp/2.0_dp))

!noise=noiseligo

!*******
!LIGO
!*******

aa0=47.8466_dp
aa1=-92.1896_dp
aa2=35.9273_dp
aa3=-7.61447_dp
aa4=0.916742_dp
aa5=-0.0588089_dp
aa6=0.00156345_dp

f0=log(f)

noiseligo=exp(aa0 + aa1*f0 + aa2*f0**2 + aa3*f0**3 &
          + aa4*f0**4 + aa5*f0**5 + aa6*f0**6)

noise=noiseligo**2

!call wave(f,mc,eta,chi,dl,tc,phic,ad1,add1,f3)

!integ=integ+amp**2/noise*w(n)
  
do i=1,imax
  do j=1,imax
     fm(i,j)=fm(i,j) &
             +real((dad1(i)*conjg(dad1(j)))/noise) &
             *4.0_dp*w(n)*ndet
  end do
end do
  
! sn=sn+ad1**2/noiseligo*4.0_dp*w(n)

!write (1,'(2E22.14)') real(dad1(4)),real(ad1)

end do

!print *, sqrt(sn)
!print *, dad1(1)
!print *, dad1(2)
!print *, dad1(3)
!print *, dad1(4)
!print *, dad1(5)
!print *, dad1(6)
!print *, dad1(7)
!print *, ad1

!do i=1,imax
!  fm(i,5)=0.0_dp
!  fm(5,i)=0.0_dp
!  fm(i,6)=0.0_dp
!  fm(6,i)=0.0_dp
!  fm(i,3)=0.0_dp
!  fm(3,i)=0.0_dp
!  fm(i,4)=0.0_dp
!  fm(4,i)=0.0_dp
!end do

!fm(5,5)=1.0_dp
!fm(6,6)=1.0_dp
!fm(3,3)=1.0_dp
!fm(4,4)=1.0_dp

!fm=fm*SNR_des**2.0_dp/fm(4,4)

fmi=fm

!print *, fm(4,4)

call gaussj(fmi,imax,imax,b,0,0)


!sn=2.0_dp*amp/dsub*sqrt(integ)

!print *, sn

end subroutine fish


subroutine wave(f,mc,eta,chi,const,tc,phic,betag,ad1,add1,f3)

use nrtype_qp
implicit none

real(dp), intent(in) :: f,mc,eta,chi,const,tc,phic,betag
real(dp), intent(inout) :: f3
complex(dpc), intent(inout) :: ad1,add1
integer :: i,j,n
real(dp) :: alpha2,alpha3,eps1,eps2,vf1,vf2,f1,f2,sigma,wm,wr,vf,psi,phase
complex(dpc) :: wi=(0.0_dp, 1.0_dp)
real(dp), dimension(4,7) :: y
real(dp), dimension(5,7) :: z
real(dp), dimension(5) :: phi

mtotal=eta**(-3.0_dp/5.0_dp)*mc


alpha2=-323.0_dp/224.0_dp+451.0_dp*eta/168.0_dp
alpha3=(27.0_dp/8.0_dp-11.0_dp*eta/6.0_dp)*chi

eps1=1.4547_dp*chi-1.8897_dp
eps2=-1.8153_dp*chi+1.6557_dp

!LIGO
!const=(mtotal*eta**(3.0_dp/5.0_dp))**(5.0_dp/6.0_dp) &
!       /(sqrt(30.0_dp)*pi_d**(2.0_dp/3.0_dp))/dl
!DPF
!const=(mtotal*eta**(3.0_dp/5.0_dp))**(5.0_dp/6.0_dp) &
!       /(pi_d**(2.0_dp/3.0_dp))/dl*sqrt(5.0_dp/24.0_dp)

!DPF 1arm
!const=mc**(5.0_dp/6.0_dp)/pi_d**(2.0_dp/3.0_dp)*sqrt(10)/15



y(1,1)=1.0_dp-4.455_dp*(1.0_dp-chi)**(0.217_dp)+3.521_dp*(1.0_dp-chi)**(0.26_dp)
y(2,1)=(1.0_dp-0.63_dp*(1.0_dp-chi)**(0.3_dp))/2.0_dp
y(3,1)=(1.0_dp-0.63_dp*(1.0_dp-chi)**(0.3_dp))*(1.0_dp-chi)**(0.45_dp)/4.0_dp
y(4,1)=0.3236_dp+0.04894_dp*chi+0.01346_dp*chi**(2.0_dp)

y(1,2)=0.6437_dp
y(1,3)=0.827_dp
y(1,4)=-0.2706_dp
y(1,5)=-0.05822_dp
y(1,6)=-3.935_dp
y(1,7)=-7.092_dp

y(2,2)=0.1469_dp
y(2,3)=-0.1228_dp
y(2,4)=-0.02609_dp
y(2,5)=-0.0249_dp
y(2,6)=0.1701_dp
y(2,7)=2.325_dp

y(3,2)=-0.4098_dp
y(3,3)=-0.03523_dp
y(3,4)=0.1008_dp
y(3,5)=1.829_dp
y(3,6)=-0.02017_dp
y(3,7)=-2.87_dp

y(4,2)=-0.1331_dp
y(4,3)=-0.08172_dp
y(4,4)=0.1451_dp
y(4,5)=-0.2714_dp
y(4,6)=0.1279_dp
y(4,7)=4.922_dp


z(1,1)=3715.0_dp/756.0_dp
z(2,1)=-16.0_dp*pi_d+113.0_dp*chi/3.0_dp
z(3,1)=15293365.0_dp/508032.0_dp-405.0_dp*chi**(2.0_dp)/8.0_dp
z(4,1)=0.0_dp
z(5,1)=0.0_dp

z(1,2)=-920.9_dp
z(1,3)=492.1_dp
z(1,4)=135.0_dp
z(1,5)=6742.0_dp
z(1,6)=-1053.0_dp
z(1,7)=-1.34_dp*10.0_dp**4

z(2,2)=1.702_dp*10.0_dp**4
z(2,3)=-9566.0_dp
z(2,4)=-2182.0_dp
z(2,5)=-1.214_dp*10.0_dp**5
z(2,6)=2.075_dp*10.0_dp**4
z(2,7)=2.386_dp*10.0_dp**5

z(3,2)=-1.254_dp*10.0_dp**5
z(3,3)=7.507_dp*10.0_dp**4
z(3,4)=1.338_dp*10.0_dp**4
z(3,5)=8.735_dp*10.0_dp**5
z(3,6)=-1.657_dp*10.0_dp**5
z(3,7)=-1.694_dp*10.0_dp**6

z(4,2)=-8.898_dp*10.0_dp**5
z(4,3)=6.31_dp*10.0_dp**5
z(4,4)=5.068_dp*10.0_dp**4
z(4,5)=5.981_dp*10.0_dp**6
z(4,6)=-1.415_dp*10.0_dp**6
z(4,7)=-1.128_dp*10.0_dp**7

z(5,2)=8.696_dp*10.0_dp**5
z(5,3)=-6.71_dp*10.0_dp**5
z(5,4)=-3.008_dp*10.0_dp**4
z(5,5)=-5.838_dp*10.0_dp**6
z(5,6)=1.514_dp*10.0_dp**6
z(5,7)=1.089_dp*10.0_dp**7


f1=(y(1,1)+y(1,2)*eta+y(1,3)*eta*chi+y(1,4)*eta*chi**2 &
    +y(1,5)*eta**2+y(1,6)*eta**2*chi+y(1,7)*eta**3)/(pi_d*mtotal)
f2=(y(2,1)+y(2,2)*eta+y(2,3)*eta*chi+y(2,4)*eta*chi**2 &
    +y(2,5)*eta**2+y(2,6)*eta**2*chi+y(2,7)*eta**3)/(pi_d*mtotal)
sigma=(y(3,1)+y(3,2)*eta+y(3,3)*eta*chi+y(3,4)*eta*chi**2 &
    +y(3,5)*eta**2+y(3,6)*eta**2*chi+y(3,7)*eta**3)/(pi_d*mtotal)
f3=(y(4,1)+y(4,2)*eta+y(4,3)*eta*chi+y(4,4)*eta*chi**2 &
    +y(4,5)*eta**2+y(4,6)*eta**2*chi+y(4,7)*eta**3)/(pi_d*mtotal)


vf1=(pi_d*mtotal*f1)**(1.0_dp/3.0_dp)
vf2=(pi_d*mtotal*f2)**(1.0_dp/3.0_dp)

vf=(pi_d*mtotal*f)**(1.0_dp/3.0_dp)

wm=(1.0_dp+alpha2*vf1**2+alpha3*vf1**3)/(1.0_dp+eps1*vf1+eps2*vf1**2)

wr=pi_d*sigma/2.0_dp*(f2/f1)**(-2.0_dp/3.0_dp)*wm*(1.0_dp+eps1*vf2+eps2*vf2**2)


if (f < f1) then
        amp = const*f**(-7.0_dp/6.0_dp)*(1.0_dp+alpha2*vf**2+alpha3*vf**3)
        
  else if (f < f2) then
        amp = const*f1**(-7.0_dp/6.0_dp)*wm*(f/f1)**(-2.0_dp/3.0_dp) &
             *(1.0_dp+eps1*vf+eps2*vf**2)
       
  else if (f < f3) then 
        amp = const*f1**(-7.0_dp/6.0_dp)*wr*sigma/(2.0_dp*pi_d) &
             /((f-f2)**2+sigma**2/4.0_dp)
         
  else  
        amp = 0.0_dp
  
end if

do i=1,5

phi(i)=z(i,1)+(z(i,2)+z(i,3)*chi+z(i,4)*chi**2)*eta+(z(i,5)+z(i,6)*chi)*eta**2+z(i,7)*eta**3

end do


if (k==2) then

psi=2.0_dp*pi_d*f*tc+phic+3.0_dp/(128.0_dp*eta*vf**5) &
*(1.0_dp+phi(1)*vf**2+phi(2)*vf**3+phi(3)*vf**4+phi(4)*vf**6+phi(5)*vf**7) &
+betag*log(pi_d*mtotal*eta**(3._dp/5._dp)*f)


else

psi=2.0_dp*pi_d*f*tc+phic+3.0_dp/(128.0_dp*eta*vf**5) &
    *(1.0_dp+phi(1)*vf**2+phi(2)*vf**3+phi(3)*vf**4+phi(4)*vf**6+phi(5)*vf**7) &
    +betag*(pi_d*mtotal*eta**(3._dp/5._dp)*f)**(b)

end if

ad1=amp*exp(wi*psi)
!ad1=const/f**(7.0_dp/6.0_dp)*exp(wi*psi)

!phase=2.0_dp*pi_d*f*tc-phic-pi_d/4.0_dp &
!            +3.0_dp/(128.0_dp*eta*vf**5) &
!            *(1.0_dp &
!            +20.0_dp/9.0_dp*(743.0_dp/336.0_dp+11.0_dp/4.0_dp*eta)*vf*2 &
!            +(-16.0_dp*pi_d)*vf**(3) & 
!            +(15293365.0_dp/508032.0_dp &
!            +27145.0_dp/504.0_dp*eta+3085.0_dp/72.0_dp*eta**2)*vf**4)

!ad1=const/f**(7.0_dp/6.0_dp)*exp(wi*phase)

end subroutine wave




SUBROUTINE gauleg(x1,x2,x,w,n)

!   INTEGER, PARAMETER :: n
INTEGER :: i,j,m,n
REAL(DP) :: x1,x2,EPS=3.q-14,p1,p2,p3,pp,xl,xm,z,z1
REAL(DP), DIMENSION(n) :: x,w
! DOUBLE PRECISION EPS
! PARAMETER (EPS=3.q-14)

!  DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1

m=(n+1)/2
xm=0.5q0*(x2+x1)
xl=0.5q0*(x2-x1)
do 12 i=1,m
z=cos(3.141592654q0*(i-.25q0)/(n+.5q0))
1       continue
p1=1.q0
p2=0.q0
do 11 j=1,n
p3=p2
p2=p1
p1=((2.q0*j-1.q0)*z*p2-(j-1.q0)*p3)/j
11        continue
pp=n*(z*p1-p2)/(z*z-1.q0)
z1=z
z=z1-p1/pp
if(abs(z-z1).gt.EPS)goto 1
x(i)=xm-xl*z
x(n+1-i)=xm+xl*z
w(i)=2.q0*xl/((1.q0-z*z)*pp*pp)
w(n+1-i)=w(i)
12    continue
return
END SUBROUTINE



SUBROUTINE gaussj(a,n,np,b,m,mp)

use nrtype_qp
implicit none

INTEGER :: m,mp,n,np
INTEGER, PARAMETER :: NMAX=50
REAL(DP), DIMENSION(np,np) :: a
REAL(DP), DIMENSION(np,mp) :: b
!  PARAMETER (NMAX=50)
INTEGER :: i,icol,irow,j,k,l,ll
INTEGER, DIMENSION(NMAX) :: indxc,indxr,ipiv
REAL(DP) :: big,dum,pivinv
do 11 j=1,n
ipiv(j)=0
11    continue
do 22 i=1,n
big=0.q0
do 13 j=1,n
if(ipiv(j).ne.1)then
do 12 k=1,n
if (ipiv(k).eq.0) then
if (abs(a(j,k)).ge.big)then
big=abs(a(j,k))
irow=j
icol=k
endif
else if (ipiv(k).gt.1) then
pause 'singular matrix in gaussj'
endif
12          continue
endif
13      continue
ipiv(icol)=ipiv(icol)+1
if (irow.ne.icol) then
do 14 l=1,n
dum=a(irow,l)
a(irow,l)=a(icol,l)
a(icol,l)=dum
14        continue
do 15 l=1,m
dum=b(irow,l)
b(irow,l)=b(icol,l)
b(icol,l)=dum
15        continue
endif
indxr(i)=irow
indxc(i)=icol
if (a(icol,icol).eq.0.q0) pause 'singular matrix in gaussj'
pivinv=1.q0/a(icol,icol)
a(icol,icol)=1.q0
do 16 l=1,n
a(icol,l)=a(icol,l)*pivinv
16      continue
do 17 l=1,m
b(icol,l)=b(icol,l)*pivinv
17      continue
do 21 ll=1,n
if(ll.ne.icol)then
dum=a(ll,icol)
a(ll,icol)=0.q0
do 18 l=1,n
a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
do 19 l=1,m
b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
endif
21      continue
22    continue
do 24 l=n,1,-1
if(indxr(l).ne.indxc(l))then
do 23 k=1,n
dum=a(k,indxr(l))
a(k,indxr(l))=a(k,indxc(l))
a(k,indxc(l))=dum
23        continue
endif
24    continue
return

END subroutine gaussj

end program bd_lisa_noangle
