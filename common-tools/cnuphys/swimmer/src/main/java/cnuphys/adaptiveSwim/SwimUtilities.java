package cnuphys.adaptiveSwim;

import cnuphys.rk4.IDerivative;

public class SwimUtilities {

	

	/**
	 * Take a single step using basic fourth order RK
	 * 
	 * @param s     the independent variable
	 * @param u     the current state vector
	 * @param du    the current derivatives
	 * @param h     the step size
	 * @param deriv can compute the rhs of the diffy q
	 * @param uf    the state vector after the step
	 */
	public static void singleRK4Step(double s, double[] u, double[] du, double h, IDerivative deriv, double[] uf) {

		int nDim = u.length;

		// note that du (input) is k1
		double k1[] = du; // the current derivatives

		double k2[] = new double[nDim];
		double k3[] = new double[nDim];
		double k4[] = new double[nDim];
		double utemp[] = new double[nDim];

		double hh = h * 0.5; // half step
		double h6 = h / 6.0;

		// advance t to mid point
		double sMid = s + hh;

		// first step: initial derivs to midpoint
		for (int i = 0; i < nDim; i++) {
			utemp[i] = u[i] + hh * k1[i];
		}
		deriv.derivative(sMid, utemp, k2);

		// 2nd step (like 1st, but use midpoint just computed derivs dyt)
		for (int i = 0; i < nDim; i++) {
			utemp[i] = u[i] + hh * k2[i];
		}
		deriv.derivative(sMid, utemp, k3);

		// third (full) step
		for (int i = 0; i < nDim; i++) {
			utemp[i] = u[i] + h * k3[i];
		}
		deriv.derivative(s + h, utemp, k4);

		for (int i = 0; i < nDim; i++) {
			uf[i] = u[i] + h6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
		}

	}
}
