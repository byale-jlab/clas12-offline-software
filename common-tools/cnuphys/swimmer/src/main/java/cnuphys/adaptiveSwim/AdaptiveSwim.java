package cnuphys.adaptiveSwim;

import cnuphys.magfield.FastMath;
import cnuphys.magfield.FieldProbe;
import cnuphys.rk4.IDerivative;
import cnuphys.rk4.IRkListener;
import cnuphys.rk4.RungeKuttaException;
import cnuphys.swim.DefaultDerivative;
import cnuphys.swim.SwimResult;

public class AdaptiveSwim {
	
	// Min momentum to swim in GeV/c
	public static final double MINMOMENTUM = 5e-05;

	//tolerance when swimmimg to a max path length
	public static final double SMAX_TOLERANCE = 1.0e-4;  //meters

	
	private static final double _minStepSize = 1.0e-5;
	private static final double _maxStepSize = 0.4;


	public static void swimRho(int charge, double xo, double yo, double zo, double momentum, double theta, double phi,
			final double targetRho, double accuracy, double s0, double sf, double h, double eps, SwimResult result)
			throws RungeKuttaException {
		
		// set u to the starting state vector
		double thetaRad = Math.toRadians(theta);
		double phiRad = Math.toRadians(phi);
		double sinTheta = Math.sin(thetaRad);

		double px = sinTheta*Math.cos(phiRad); //px/p
		double py = sinTheta*Math.sin(phiRad); //py/p
		double pz = Math.cos(thetaRad); //pz/p
		
		double uf[] = result.getUf();
		uf[0] = xo;
		uf[1] = yo;
		uf[2] = xo;
		uf[3] = px;
		uf[4] = py;
		uf[5] = pz;

		if (momentum < MINMOMENTUM) {
			System.err.println("Skipping low momentum fixed rho swim (A)");
			result.setNStep(0);
			result.setFinalS(0);
			result.setStatus(-2);
			return;
		}
		
		//cutoff value of s with tolerance 
		double sCutoff = sf - SMAX_TOLERANCE;
		
		double del = Double.POSITIVE_INFINITY;
		int maxtry = 11;
		int count = 0;
		double sCurrent = 0;
		int ns = 0;

		//swimmer with current field
		FieldProbe probe = FieldProbe.factory();
		DefaultDerivative deriv = new DefaultDerivative(charge, momentum, probe);

		AdaptiveRhoStopper stopper = new AdaptiveRhoStopper(uf.length, targetRho, accuracy);

		//use a half-step advancer
		HalfStepAdvance advance = new HalfStepAdvance();

		while ((count < maxtry) && (del > accuracy)) {
			
			uf = result.getUf();
			if (count > 0) {
				px = uf[3];
				py = uf[4];
				pz = uf[5];
				theta = FastMath.acos2Deg(pz);
				phi = FastMath.atan2Deg(py, px);
			}
			
			double rhoCurr = FastMath.hypot(uf[0], uf[1]);
			
			stopper.set(sCurrent, sf, rhoCurr);
			
			if ((sCurrent + h) > sf) {
				h = (sf-sCurrent)/4;
				if (h < 0) {
					break;
				}
			}

			ns += driver(uf, 0, sf-sCurrent, h, deriv, stopper, null, advance, eps, uf);
			
			System.arraycopy(stopper.getFinalU(), 0, result.getUf(), 0, result.getUf().length);
			
			sCurrent = stopper.getFinalS();
			
			double rholast = FastMath.hypot(result.getUf()[0], result.getUf()[1]);
			del = Math.abs(rholast - targetRho);
			
			if ((sCurrent) > sCutoff) {
//				System.out.println(" s final   " + (finalPathLength + stepSize) +  "  rhoLast = " + rholast);
				break;
			}
			
			count++;
			h = Math.min(h, (sf-sCurrent)/4);
//			stepSize /= 2;


		}
		
		
		result.setFinalS(sCurrent);
		result.setUf(uf);
		result.setNStep(ns);

		if (del < accuracy) {
			result.setStatus(0);
		} else {
			result.setStatus(-1);
		}
	}
	
	
	/**
	 * Driver that uses the RungeKutta advance with an adaptive step size
	 * 
	 * This version uses an IRk4Listener to notify the listener that the next step
	 * has been advanced.
	 * 
	 * A very typical case is a 2nd order ODE converted to a 1st order where the
	 * dependent variables are x, y, z, px/p, py/p, pz/p and the independent
	 * variable is time.
	 * 
	 * @param uo           initial values. Probably something like (xo, yo, zo, vxo,
	 *                     vyo, vzo).
	 * @param s0           the initial value of the independent variable, e.g.,
	 *                     time.
	 * @param sf           the maximum value of the independent variable.
	 * @param h            the step size
	 * @param deriv        the derivative computer (interface). This is where the
	 *                     problem specificity resides.
	 * @param stopper      if not <code>null</code> will be used to exit the
	 *                     integration early because some condition has been
	 *                     reached.
	 * @param listener     listens for each step
	 * @param advancer     takes the next step
	 * @param eps          sort of relative tolerance (e.g., 1.0e-6)
	 * @param uf           will hold final state vector
	 *                     
	 * @return the number of steps used.
	 * @throw(new RungeKuttaException("Step size too small in Runge Kutta driver"
	 *            ));
	 */
	public static int driver(double uo[], double s0, double sf, double h, IDerivative deriv, IAdaptiveStopper stopper,
			IRkListener listener, IAdaptiveAdvance advancer, double eps, double uf[]) throws RungeKuttaException {

		// the dimensionality of the problem. E.., 6 if (x, y, z, vx, vy, vz)
		int nDim = uo.length;
		
		//running val of independent variable
		double s = s0;
		
		// ut is the running value of the state vector,
		// typically [x, y, z, px/p, py/p, pz/p]
		double ut[] = new double[nDim];
		
		//du is for derivatives
		double du[] = new double[nDim];

		//init the final value at the start value
		System.arraycopy(uo, 0, uf, 0, nDim);
		
		int nstep = 0;
		
		//we are done when this sign changes (or if the stopper
		//stops us first)
		int startSign = sign(sf-s0);
		
		//no integration range??
		if (startSign == 0) {
			return 0;
		}
		
		AdaptiveStepResult result = new AdaptiveStepResult();
		
		while (sign(sf-s) == startSign) {
			
			System.arraycopy(uf, 0, ut, 0, nDim);

			//compute derivs at current step
			deriv.derivative(s, ut, du);
			advancer.advance(s, ut, du, h, deriv, uf, eps, result);
			
			h = result.getHNew();
			s = result.getSNew();
			
			if (Math.abs(h) < _minStepSize) {
				h = _minStepSize;
//				throw (new RungeKuttaException("Step size too small in AdaptiveSwim driver"));
			}
			else {
				nstep++;
				
				// someone listening?
				if (listener != null) {
					listener.nextStep(s, uf, h);
				}
				
				if (stopper.stopIntegration(s, uf)) {
					return nstep;
				}


			}

		}
		
		return nstep;
		

	}


	//sign function
	private static int sign(double v) {
		if (v > 0) {
			return 1;
		}
		else if (v < 0) {
			return -1;
		}
		return 0;
	}
	
   		
}
