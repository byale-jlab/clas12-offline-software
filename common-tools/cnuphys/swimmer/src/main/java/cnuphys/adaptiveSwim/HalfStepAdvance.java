package cnuphys.adaptiveSwim;

import cnuphys.rk4.IDerivative;

public class HalfStepAdvance implements IAdaptiveAdvance {

	private static final double _safety = 0.95;
	private static final double _pgrow = -0.20;
	private static final double _pshrink = -0.25;
	private static final double _errControl = 1.0e-4;
	private static final double _correctFifth = 1. / 15;
	private static final double _tiny = 1.0e-30;

	@Override
	public void advance(double s, double[] u, double[] du, double h, IDerivative deriv, double[] uf, double eps,
			AdaptiveStepResult result) {

		boolean done = false;

		int ndim = u.length;
		double usave[] = new double[ndim];
		double dusave[] = new double[ndim];
		double utemp[] = new double[ndim];
		double uscale[] = new double[ndim];

		while (!done) {
			//save in case our stepsize was not accepted (too big)
			System.arraycopy(u, 0, usave, 0, ndim);
			System.arraycopy(du, 0, dusave, 0, ndim);

			// almost relative error, but with safety when values of u are small
			for (int i = 0; i < ndim; i++) {
				uscale[i] = Math.abs(usave[i]) + Math.abs(h * dusave[i]) + _tiny;
			}

			// advance two half steps after which uf will hold the value of 
			// which, if our steps size is acceptable, this will be our result
			double h2 = h / 2;
			double smid = s + h2;
			SwimUtilities.singleRK4Step(s, usave, dusave, h2, deriv, utemp);
			deriv.derivative(smid, utemp, du);
			SwimUtilities.singleRK4Step(smid, utemp, du, h2, deriv, uf);

			// take the full step
			SwimUtilities.singleRK4Step(s, usave, dusave, h, deriv, utemp);

			// compute the maximum error
			double errMax = 0;
			for (int i = 0; i < ndim; i++) {
				utemp[i] = uf[i] - utemp[i];
				errMax = Math.max(errMax, Math.abs(utemp[i] / uscale[i]));
			}

			// scale based on tolerance in eps
			errMax = errMax / eps;

			if (errMax > 1) {
				//get smaller h, then try again since done = false
				
				double shrinkFact = _safety * Math.pow(errMax, _pshrink);
				
				//no more than a factor of 4
	//			shrinkFact = Math.max(shrinkFact, 0.25);
				h = h * shrinkFact;
				

			} else { // can grow
				double hnew;
				if (errMax > _errControl) {
					hnew = _safety * h * Math.pow(errMax, _pgrow);
				} else {
					hnew = 4 * h;
				}
				
				result.setHNew(hnew);
				result.setHUsed(h);
				result.setSNew(s + h);
				done = true;
			}

		} // !done

		// mop up 5th order truncation error
		//so result is actuallky 5th order
		for (int i = 0; i < ndim; i++) {
			uf[i] = uf[i] + utemp[i] * _correctFifth;
		}
	} // end advance
}
