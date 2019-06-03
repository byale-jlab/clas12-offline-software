package cnuphys.swimtest;

import java.util.Random;

import cnuphys.adaptiveSwim.AdaptiveSwim;
import cnuphys.magfield.MagneticFields;
import cnuphys.magfield.MagneticFields.FieldType;
import cnuphys.rk4.RungeKuttaException;
import cnuphys.swim.SwimResult;
import cnuphys.swim.Swimmer;

public class RhoTest {

	/** Test swimming to a fixed rho (cylinder) */
	public static void rhoTest() {
		
		long seed = 9459363;
		System.out.println("TEST swimming to a fixed rho (cylinder)");
		MagneticFields.getInstance().setActiveField(FieldType.SOLENOID);
		
		double xo = 0; // m
		double yo = 0; // m
		double zo = 0; // m
		double stepsizeAdaptive = 0.01; //starting
		double stepsizeUniform = 5e-04;  //m
		
		double maxPathLength = 3; //m
		double accuracy = 5e-3; //m
		double fixedRho = 0.26;  //m 
		
		int num = 10000;
//		num = 10000;
		int n0= 0;
		
		SwimResult uniform = new SwimResult(6);
		SwimResult adaptive = new SwimResult(6);
		SwimResult newadaptive = new SwimResult(6);
		
		//generate some random initial conditions
		Random rand = new Random(seed);
		
		int charge[] = new int[num];
		double p[] = new double[num];
		double theta[] = new double[num];
		double phi[] = new double[num];
		int adaptStatus[] = new int[num];
		
		for (int i = 0; i < num; i++) {
			charge[i] = (rand.nextDouble() < 0.5) ? -1 : 1;
//			p[i] = 1 + 8*rand.nextDouble();
//			theta[i] = 45 + 25*rand.nextDouble();
//			phi[i] = 360*rand.nextDouble();
			p[i] = 0.25 + 0.75*rand.nextDouble();
			theta[i] = 40 + 30*rand.nextDouble();
			phi[i] = 360*rand.nextDouble();
		}
		
		
		long time;
        double rhof;
		double sum;
		double delMax;
		Swimmer swimmer = new Swimmer();
		int badStatusCount;
		
		long nStepTotal = 0;
		
		// adaptive step
		try {
			
			sum = 0;
			badStatusCount = 0;
			delMax = Double.NEGATIVE_INFINITY;
			time = System.currentTimeMillis();

			
			for (int i = n0; i < num; i++) {
				swimmer.swimRho(charge[i], xo, yo, zo, p[i], theta[i], phi[i], fixedRho, accuracy, maxPathLength,
						stepsizeAdaptive, Swimmer.CLAS_Tolerance, adaptive);

				rhof = Math.hypot(adaptive.getUf()[0], adaptive.getUf()[1]);
				double dd = Math.abs(fixedRho - rhof);
				delMax = Math.max(delMax, dd);
				sum += dd;
				
				adaptStatus[i] = adaptive.getStatus();
				nStepTotal += adaptive.getNStep();
				
				if (adaptive.getStatus() != 0) {
					badStatusCount += 1;
	//				num = i+1;
				}
			}
			
			time = System.currentTimeMillis() - time;
			SwimTest.printSummary("Fixed Rho,  Adaptive step size", adaptive.getNStep(), p[num - 1],
					adaptive.getUf(), null);
			System.out.println(String.format("Adaptive time: %-7.3f   avg delta = %-9.5f  max delta = %-9.5f  badStatCnt = %d", ((double)time)/1000., sum/num, delMax, badStatusCount));
			System.out.println("Adaptive Avg NS = " +  (int)(((double)nStepTotal)/num));
			System.out.println("Adaptive Path length = " + adaptive.getFinalS() + " m\n\n");

		} catch (RungeKuttaException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		// NEW adaptive step
		try {
			
			sum = 0;
			badStatusCount = 0;
			delMax = Double.NEGATIVE_INFINITY;
			time = System.currentTimeMillis();
			double eps = 1.0e-6;
			nStepTotal = 0;

			for (int i = n0; i < num; i++) {

				AdaptiveSwim.swimRho(charge[i], xo, yo, zo, p[i], theta[i], phi[i], fixedRho, accuracy, 0, maxPathLength, stepsizeAdaptive, eps, newadaptive);

				rhof = Math.hypot(newadaptive.getUf()[0], newadaptive.getUf()[1]);
				double dd = Math.abs(fixedRho - rhof);
				delMax = Math.max(delMax, dd);
				sum += dd;
				
				if (newadaptive.getStatus() != adaptStatus[i]) {
					System.out.println("Adaptive v. NEW Adaptive Status differs for i = " + i + "     adaptiveStat = " + adaptStatus[i]
							+ "    NEW adaptive status = " + newadaptive.getStatus());
				}
				
				nStepTotal += newadaptive.getNStep();
				if (newadaptive.getStatus() != 0) {
					badStatusCount += 1;
	//				num = i+1;
				}
			}
			
			time = System.currentTimeMillis() - time;
			SwimTest.printSummary("Fixed Rho,  NEW Adaptive step size", newadaptive.getNStep(), p[num - 1],
					newadaptive.getUf(), null);
			System.out.println(String.format("NEW Adaptive time: %-7.3f   avg delta = %-9.5f  max delta = %-9.5f  badStatCnt = %d", ((double)time)/1000., sum/num, delMax, badStatusCount));
			System.out.println("NEW Adaptive Avg NS = " +  (int)(((double)nStepTotal)/num));
			System.out.println("NEW Adaptive Path length = " + newadaptive.getFinalS() + " m\n\n");

		} catch (RungeKuttaException e) {
			e.printStackTrace();
		}


		// uniform step
		time = System.currentTimeMillis();

		sum = 0;
		badStatusCount = 0;
		nStepTotal = 0;
		delMax = Double.NEGATIVE_INFINITY;
		for (int i = n0; i < num; i++) {
			swimmer.swimRho(charge[i], xo, yo, zo, p[i], theta[i], phi[i], fixedRho, accuracy, maxPathLength, stepsizeUniform, uniform);			
			rhof = Math.hypot(uniform.getUf()[0], uniform.getUf()[1]);
			double dd = Math.abs(fixedRho - rhof);
			delMax = Math.max(delMax, dd);
			sum += dd;
			
			if (uniform.getStatus() != adaptStatus[i]) {
				System.out.println("Status differs for i = " + i + "     adaptiveStat = " + adaptStatus[i] + "    uniform status = " + uniform.getStatus());
			}

			nStepTotal += uniform.getNStep();
			if (uniform.getStatus() != 0) {
				badStatusCount += 1;
			}
		}
		time = System.currentTimeMillis() - time;
		SwimTest.printSummary("Fixed Rho,  Uniform step size", uniform.getNStep(), p[num - 1],
				uniform.getUf(), null);
		System.out.println(String.format("Uniform time: %-7.3f   avg delta = %-9.5f  max delta = %-9.5f  badStatCnt = %d", ((double)time)/1000., sum/num, delMax, badStatusCount));
		System.out.println("Uniform Avg NS = " +  (int)(((double)nStepTotal)/num));
		System.out.println("Uniform Path length = " + uniform.getFinalS() + " m\n\n");
		
		
		
		
	}

}
