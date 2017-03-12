#! /usr/bin/env python
# -*- coding: utf-8 -*-


from libv import *
import math
import matplotlib.pyplot as plt
import numpy as np
import unittest
import scipy


np.set_printoptions(precision=6)
np.seterr(divide='ignore', invalid='ignore')


class vawt_test(unittest.TestCase):
    def test_airfoil_extrapolate(self):
       print "test_airfoil_extrapolate"
       domain = libv()
       f1 = domain.airfoil("../foils/NACA_0015_extrap.dat",name="foil1")        
       f2 = domain.airfoil("../foils/NACA_0015.dat",name="foil2")
       f3 = domain.airfoil("../foils/NACA_0015.dat",name="foil2_extrap")
       f3.extrapolate(1.2)    
       int_cl1 = scipy.integrate.simps(f1.cl, f1.alpha)
       int_cl3 = scipy.integrate.simps(f3.cl, f3.alpha)        
       self.assertLessEqual(int_cl1-int_cl3, 1E-7, "integrared areas for foil 1 and foil 3 are not identical")
       
       plt.figure(1)        
       plt.subplot(2,1,1)
       plt.plot(f1.alpha, f1.cl,'r-',lw=6,alpha=0.25,label='Foil 1')
       plt.plot(f2.alpha, f2.cl,'g--',lw=4,alpha=0.5,label='Foil 2')
       plt.plot(f3.alpha, f3.cl,'k-',lw=1,label='Extrapolated')
       # plt.xlabel('Angle of Attack [rad]',fontsize=14)
       plt.ylabel('$C_L$',fontsize=14)
       plt.grid(True)
       plt.legend(fontsize=10,loc=2)
       
       plt.subplot(2,1,2)
       plt.plot(f1.alpha, f1.cd,'r-',lw=6,alpha=0.25,label='Foil 1')
       plt.plot(f2.alpha, f2.cd,'g--',lw=4,alpha=0.5,label='Foil 2')
       plt.plot(f3.alpha, f3.cd,'k-',lw=1,label='Extrapolated')
       plt.xlabel('Angle of Attack [rad]',fontsize=14)
       plt.ylabel('$C_D$',fontsize=14)
       plt.grid(True)
       plt.legend(fontsize=10,loc=2)

       # plt.savefig('test_airfoil_extrapolate.png', dpi=150)
       plt.show()
       
       
    def test_airfoil_blend(self):
       print "test_airfoil_blend"
       domain = libv()
       f1 = domain.airfoil("../foils/NACA_0015_extrap.dat",name="foil2")
       f2 = domain.airfoil("../foils/NACA64_A17.dat",name="foil2212121")
       f3 = domain.airfoil("../foils/NACA_0015.dat",name="foil1")
       blend = domain.airfoil_blend(0.5, f1, f2, name="BLENDED")
   
       plt.figure(2)
       plt.subplot(2,1,1)
       plt.plot(f1.alpha, f1.cl,'r-',lw=4,alpha=0.5,label='Foil A Extrapolated')
       plt.plot(f2.alpha, f2.cl,'k-',lw=1,label='Foil B')
       plt.plot(f3.alpha, f3.cl,'g--',lw=3,label='Foil A')
       plt.plot(blend.alpha, blend.cl,'b',lw=4,alpha=0.5, label='Blended (50%)')
       plt.ylabel('$C_L$',fontsize=14)
       plt.grid(True)
       plt.legend(fontsize=10,loc=2)
       
       plt.subplot(2,1,2)
       plt.plot(f1.alpha, f1.cd,'r-',lw=4,alpha=0.5,label='Foil A Extrapolated')
       plt.plot(f2.alpha, f2.cd,'k-',lw=1,label='Foil B')
       plt.plot(f3.alpha, f3.cd,'g--',lw=3,label='Foil A')
       plt.plot(blend.alpha, blend.cd,'b',lw=4,alpha=0.5, label='Blended (50%)')
       plt.xlabel('Angle of Attack [rad]',fontsize=14)
       plt.ylabel('$C_D$',fontsize=14)
       plt.grid(True)
       plt.legend(fontsize=10,loc=2)

       # plt.savefig('test_airfoil_blend.png', dpi=150)
       plt.show()
   
       
    def test_dynamic_stall_1(self):
       print "test_dynamic_stall_1"
       domain = libv()
       H = 5.0
       hubHt = 6.0
       D = 6.0
       h = 10
       R = D/2.0        
       chord_width = 0.25
       env = domain.environment()
       env.rho = 1.225 
       env.mu = 1.7894E-5 
       env.shear = 0.185 
       env.ref_h = hubHt 
   
       z = np.linspace(hubHt-H/2.0, hubHt+H/2.0, h)
       r = D/2.0*np.ones_like(z)
       chord = chord_width*np.ones_like(z)
       swept = np.zeros_like(z)
       # swept = np.linspace(-2*np.pi/12, 2*np.pi/12, h)
       chord[z > hubHt+H/2.0-1] = np.linspace(chord_width, 0.6*chord_width, len(chord[z > hubHt+H/2.0-1]))
       chord[z < hubHt-H/2.0+1] = np.linspace(0.6*chord_width, chord_width, len(chord[z < hubHt-H/2.0+1]))
       twist = np.zeros_like(r)
   
       foil = domain.airfoil("../foils/NACA_0015_extrap.dat",name="foil")
       af = [foil.af_get() for i in xrange(len(r))]      
   
       vawt = domain.vawt(3, z, 36)
       vawt.set_blade_foils(0, af, chord, r, twist, swept)
       vawt.set_blade_foils(1, af, chord, r, twist, swept)
       vawt.set_blade_foils(2, af, chord, r, twist, swept)
       vawt.dstall = True
       vawt.set_initial_theta(np.pi/2)
       vawt.commit(env)
       
       domain.render_engine()
       domain.render()
   
       vawt.omega = 127.0*math.pi/30.0  # rad/s
       tsr = np.linspace(0.5, 8, 25)
       cp = []
       ct = []
       for i in range(len(tsr)): 
           Uinf = vawt.omega*R/tsr[i]
           env.u_inf = Uinf
           vawt.solve()
           vawt.performance()
           cp.append(vawt.Cp)
           ct.append(vawt.Ct)
       domain.end()
   
       cp_original = np.array([ 0.004123,  0.012306,  0.041286,  0.067983,  0.104686,  0.155869,  0.224918,
                                0.319169,  0.388652,  0.411041,  0.398778,  0.372783,  0.343109,  0.311283,
                                0.274407,  0.231657,  0.186548,  0.139962,  0.092019,  0.042843, -0.008285,
                                -0.062901, -0.122659, -0.187782, -0.258613])                               
   
       ct_original = np.array([ 0.080434,  0.111659,  0.154766,  0.194252,  0.258171,  0.338443,  0.431183,
                                0.551965,  0.655877,  0.717531,  0.742165,  0.750798,  0.756621,  0.759269,
                                0.757969,  0.754942,  0.751264,  0.747883,  0.746862,  0.749898,  0.757027,
                                0.766309,  0.775372,  0.784147,  0.792592])
   
       error_cp = np.linalg.norm(np.array(cp)-cp_original)    
       error_ct = np.linalg.norm(np.array(ct)-ct_original)
   
       self.assertLessEqual(error_cp, 1E-5)
       self.assertLessEqual(error_ct, 1E-5)

        
    def test_plot_power(self):
        print "test_plot_power"
        domain = libv()
        H = 5.0
        hubHt = 6.0
        D = 6.0
    
        env = domain.environment()
        env.rho = 1.225 
        env.mu = 1.7894E-5 
        env.shear = 0.185 
        env.ref_h = hubHt 
     
        h = 10
        z = np.linspace(hubHt-H/2.0, hubHt+H/2.0, h)
        r = D/2.0*np.ones_like(z)
        C = 0.25
        chord = C*np.ones_like(z)
        swept = np.zeros_like(z)
        twist = np.zeros_like(r)
     
        foil = domain.airfoil("../foils/NACA_0015_extrap.dat",name="foil2")
        af = [foil.af_get() for i in xrange(len(r))]

        vawt = domain.vawt(3, z, 36)
        vawt.set_blade_foils(0, af, chord, r, twist, swept)
        vawt.set_blade_foils(1, af, chord, r, twist, swept)
        vawt.set_blade_foils(2, af, chord, r, twist, swept)
        vawt.set_initial_theta(np.pi/2)
        vawt.commit(env)    
     
        vawt.omega = 127.0*math.pi/30.0  # rad/s
        R = D/2.0

        vawt.set_initial_theta(-np.pi/2)
        vawt.dstall = True
        tsr = np.linspace(1, 8, 20)
        cp = []
        ct = []
        for i in range(len(tsr)): 
            Uinf = vawt.omega*R/tsr[i]
            env.u_inf = Uinf
            vawt.solve()
            vawt.performance()
            cp.append(vawt.Cp)
            ct.append(vawt.Ct)    
        tsr_ac = np.linspace(1, 8, 20)
        CP_ac = [0.0152, 0.0370644,0.0749492,0.1348237,0.2263019,0.3574439,0.4405602,0.4569397,0.4439087,0.4193983,0.3876984,0.3494999,0.3052150,0.2549715,0.1986731,0.1359335,0.0662743,-0.010748,-0.095550,-0.188554]
        CT_ac = [0.10916982, 0.15157266, 0.21333917, 0.30237707, 0.4326099, 0.59959357, 0.72735771, 0.80983707, 0.86073882, 0.88985662, 0.90833265, 0.92204745, 0.93285898, 0.940987, 0.94689473, 0.95125532, 0.95455376, 0.95704846, 0.95890756, 0.96024314]
        tsr_cactus = np.arange(1, 7.1, 0.5)
        CP_cactus = [0.029, 0.0708, 0.173, 0.342, 0.468, 0.515, 0.497, 0.456, 0.416, 0.356, 0.227, 0.136, -0.0591]
        tsr_cfd = np.arange(1, 8.1, 0.5)
        CP_cfd = [0.026824836, 0.063673186, 0.150848839, 0.352186123, 0.43131109, 0.439449322, 0.414917889, 0.366730595, 0.29907876, 0.214189208, 0.108863413, -0.018983, -0.16725131, -0.340271956, -0.523360716]
        CP_c = [0.115099, 0.186130, 0.287411, 0.392098, 0.450197, 0.462926, 0.454492, 0.436138, 0.412598, 0.384756, 0.352557, 0.315846, 0.274734, 0.229196, 0.179472, 0.125436, 0.066868, 0.003418, -0.065130,-0.139065]
        tsr_c = [2.000000,2.315789,2.631579, 2.947368,3.263158, 3.578947,3.894737, 4.210526,4.526316, 4.842105,5.157895, 5.473684,5.789474,6.105263,6.421053,6.736842,7.052632,7.368421,7.684211,8.000000]
        
        plt.figure(3)
        plt.subplot(2,1,1)
        plt.plot(tsr, cp, '-bs',mew=2,mfc="None",mec='b',lw=2,alpha=0.5, label='AC+Dynamic Stall')
        plt.plot(tsr_ac, CP_ac,'-ro',mew=2,mfc="None",mec='r',alpha=0.5,lw=2,label='AC')
        plt.plot(tsr_cactus, CP_cactus, '-y^',mew=2,mfc="None",mec='y',alpha=0.5,lw=2,label='CACTUS')
        plt.plot(tsr_cfd, CP_cfd, '-cp',mew=2,mfc="None",mec='c',alpha=0.5,lw=2,label='CFD')
        plt.legend()
        plt.grid(True)
        plt.legend(fontsize=10,loc=3,numpoints=1)
        plt.ylabel('$C_P$',fontsize=14)
        
        plt.subplot(2,1,2)
        plt.plot(tsr, ct, '-bs',mew=2,mfc="None",mec='b',lw=2,alpha=0.5, label='AC+Dynamic Stall')
        plt.plot(tsr_ac, CT_ac,'-ro',mew=2,mfc="None",mec='r',alpha=0.5,lw=2,label='AC')
        plt.grid(True)
        plt.legend(fontsize=10,loc=0,numpoints=1)
        plt.xlabel('Tip Speed Ratio $\lambda$',fontsize=14)
        plt.ylabel('$C_T$',fontsize=14)

        # plt.savefig('test_plot_power.png', dpi=150)
        plt.show()
        
        fx = [ 0.74109225,  0.3869074 ,  0.1951985 ,  0.19758182,  0.39278853,
               0.7479061 ,  1.19714848,  1.58890733,  1.8102441 ,  1.81097889,
               1.58709332,  1.18870357,  0.74109225,  0.3869074 ,  0.1951985 ,
               0.19758182,  0.39278853,  0.7479061 ,  1.19714848,  1.58890733,
               1.8102441 ,  1.81097889,  1.58709332,  1.18870357,  0.74109225,
               0.3869074 ,  0.1951985 ,  0.19758182,  0.39278853,  0.7479061 ,
               1.19714848,  1.58890733,  1.8102441 ,  1.81097889,  1.58709332,
               1.18870357]
        fy = [ 0.98583392,  0.6744203 ,  0.25291381, -0.20010485, -0.60947218,
               -0.91213314, -1.02756834, -0.77708787, -0.3089817 ,  0.2397277 ,
               0.74229079,  1.06102904,  0.98583392,  0.6744203 ,  0.25291381,
               -0.20010485, -0.60947218, -0.91213314, -1.02756834, -0.77708787,
               -0.3089817 ,  0.2397277 ,  0.74229079,  1.06102904,  0.98583392,
               0.6744203 ,  0.25291381, -0.20010485, -0.60947218, -0.91213314,
               -1.02756834, -0.77708787, -0.3089817 ,  0.2397277 ,  0.74229079,
               1.0610]

        error_cp = np.linalg.norm(np.array(cp)-np.array(CP_ac))      
        error_ct = np.linalg.norm(np.array(ct)-np.array(CT_ac))
        error_fx = np.linalg.norm(np.array(fx)-np.array(vawt.fx))      
        error_fy = np.linalg.norm(np.array(fy)-np.array(vawt.fy))


        self.assertAlmostEqual(error_cp, 0.387786466)
        self.assertAlmostEqual(error_ct, 0.549409777)
        self.assertAlmostEqual(error_fx, 1.173611342)
        self.assertAlmostEqual(error_fy, 1.199699893)

        
    def test_dynamic_stall_2(self):
        print "test_dynamic_stall_2"
        domain = libv()
        H = 5.0
        hubHt = 6.0
        D = 6.0
    
        env = domain.environment()
        env.rho = 1.225 
        env.mu = 1.7894E-5 
        env.shear = 0.185 
        env.ref_h = hubHt 
     
        h = 10
        z = np.linspace(hubHt-H/2.0, hubHt+H/2.0, h)
        r = D/2.0*np.ones_like(z)
        C = 0.25
        chord = C*np.ones_like(z)
        swept = np.zeros_like(z)
        twist = np.zeros_like(r)
     
        foil = domain.airfoil("../foils/NACA_0015_extrap.dat",name="foil2")
        af = [foil.af_get() for i in xrange(len(r))]
    
        vawt = domain.vawt(3, z, 36)
        vawt.set_blade_foils(0, af, chord, r, twist, swept)
        vawt.set_blade_foils(1, af, chord, r, twist, swept)
        vawt.set_blade_foils(2, af, chord, r, twist, swept)
        vawt.set_initial_theta(np.pi/2)
        vawt.commit(env)    
     
        vawt.omega = 127.0*math.pi/30.0  # rad/s
        R = D/2.0
        vawt.set_initial_theta(-np.pi/2)
    
        vawt.dstall = False
        tsr = [3]
        plt.figure(2)
        for i in range(len(tsr)): 
            Uinf = vawt.omega*R/tsr[i]
            env.u_inf = Uinf
            vawt.solve()
            vawt.performance()
            plt.subplot(2,1,1)
            plt.plot(vawt.alpha, vawt.cl, lw=2, label='Dynamic Stall = False')
            plt.legend(fontsize=10,loc=2,numpoints=1)
            plt.ylabel('$C_L$',fontsize=14)
        
            plt.subplot(2,1,2)
            plt.grid(True)
            plt.plot(vawt.alpha, vawt.cd, lw=2, label='No Dynamics Stall')
            plt.ylabel('$C_D$',fontsize=14)
            
        vawt.dstall = True
        for i in range(len(tsr)): 
            Uinf = vawt.omega*R/tsr[i]
            env.u_inf = Uinf
            vawt.solve()
            vawt.performance()
            plt.subplot(2,1,1)
            plt.grid(True)
            plt.plot(vawt.alpha, vawt.cl, 'r', lw=2, label='Unsteady Aerodynamics')
            plt.legend(fontsize=10,loc=2,numpoints=1)
            plt.ylabel('$C_L$',fontsize=14)
    
            plt.subplot(2,1,2)
            plt.grid(True)
            plt.plot(vawt.alpha, vawt.cd, 'r', lw=2, label='$\lambda=%1.1f$'%(tsr[i]))
            plt.ylabel('$C_D$',fontsize=14)
            plt.xlabel('Angle of Attack [rad]',fontsize=14)

        # plt.savefig('test_dynamic_stall_2.png', dpi=150)
        plt.show()
        
        cl = np.array([-0.14599572734162347, 0.0054064830882442445, 0.13079404730175723, 0.25676604321336643, 0.3781739439985339, 0.4941014930104249, 0.6086776120130325, 0.7020210912664334, 0.7698978340959729, 0.8082354310852236, 0.8132684441559681, 0.7816815689539283, 0.7107272299314077, 0.5988196118123705, 0.45885461613017986, 0.2950157525439928, 0.10260786063936346, -0.11211010558314041, -0.3142230165841363, -0.4366995446374996, -0.4955722433094483, -0.4988449888445762, -0.4809869489214157, -0.46549665385611666, -0.4459547890722465, -0.42465103089467215, -0.40348122789470664, -0.38530594904708526, -0.37086309273864404, -0.3610831228271349, -0.3568077259378087, -0.37651148712432003, -0.39451445501000465, -0.3943832900990085, -0.36308475589385664, -0.29095009336414807])
        cd = np.array([0.01125482494992874, 0.010711392486212862, 0.011740927969148487, 0.012972847980602473, 0.013577640146007056, 0.014269374531831396, 0.01520770642522978, 0.015342186016259519, 0.014847367200037143, 0.01371274165498847, 0.011992584967792281, 0.009867616374343333, 0.007690212996741167, 0.006007538272480396, 0.0057417581725399715, 0.006740875719066118, 0.009181537812924218, 0.012948559299535228, 0.01693626660246853, 0.015851636157636476, 0.013685423799292187, 0.011408334644922597, 0.01038914829460907, 0.01029829017407145, 0.010181391011912792, 0.010144574180770319, 0.01016625106405479, 0.010256459898807208, 0.01036238807503698, 0.010481056440769892, 0.010611182228419575, 0.01112222000331144, 0.011173973789794244, 0.010804484277070387, 0.010181784000791456, 0.009658188741964795])
        error_cl = np.array(vawt.cl)-cl
        error_cd = np.array(vawt.cd)-cd
    
        
    def test_stationary_performance(self):
        print "test_stationary_performance"
        domain = libv()
        H = 5.0
        hubHt = 6.0
        D = 6.0
    
        env = domain.environment()
        env.rho = 1.225 
        env.mu = 1.7894E-5 
        env.shear = 0.185 
        env.ref_h = hubHt 
     
        h = 10
        z = np.linspace(hubHt-H/2.0, hubHt+H/2.0, h)
        r = D/2.0*np.ones_like(z)
        C = 0.25
        chord = C*np.ones_like(z)
        swept = np.zeros_like(z)
        twist = np.zeros_like(r)
    
        foil = domain.airfoil("../foils/NACA_0015_extrap.dat",name="foil2")
        af = [foil.af_get() for i in xrange(len(r))]
    
        vawt = domain.vawt(3, z, 36)
        vawt.set_blade_foils(0, af, chord, r, twist, swept)
        vawt.set_blade_foils(1, af, chord, r, twist, swept)
        vawt.set_blade_foils(2, af, chord, r, twist, swept)
        vawt.set_initial_theta(np.pi/2)
        vawt.commit(env)    
     
        vawt.omega = 127.0*math.pi/30.0  # rad/s
        R = D/2.0
        vawt.set_initial_theta(-np.pi/2)
    
        U = np.linspace(1, 25, 5)
        theta = np.linspace(0, 2*np.pi,100)
        plt.figure(5)
        for i in range(len(U)):
            q = []
            env.u_inf = U[i]
            for j in range(len(theta)): 
                q.append(vawt.get_stationary_q(theta[j]))
            plt.plot(theta, q,label='U = %1.1f m/s'%U[i])
        plt.ylabel('Torque [N-m]',fontsize=14)
        plt.title('Torque for stationary (non-rotating) VAWT',fontsize=14)
        plt.xlabel('Rotor Offset Angle [rad]',fontsize=14)
        plt.legend(prop={'size':10})

        # plt.savefig('test_stationary_performance.png', dpi=150)
        plt.show()

        
if __name__ == '__main__':
    unittest.main()
    
