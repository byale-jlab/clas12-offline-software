package cnuphys.adaptiveSwim;

public class AdaptiveRhoStopper implements IAdaptiveStopper {

	private double _targetRho;
	private double _accuracy;
	private double _s0; //starting path length
	
	private double _totS; //total path length
	private double _prevS; //previous step path length

	private int _dim;   //dimension of our system
	private double _sMax; //max path length meters
	private double[] _uf; //final state vector

	//is the starting rho bigger or smaller than the target
	private int _startSign;
	
	/**
	 * Rho  stopper  (does check max path length)
	 * @param targetRho       stopping rho in meters
	 * @param accuracy        the accuracy in meters
	 */
	public AdaptiveRhoStopper(int dim, double targetRho, double accuracy) {
		_dim = dim;
		_targetRho = targetRho;
		_accuracy = accuracy;
		_uf = new double[_dim];
	}

		
	/**
	 * So that the stopper can be reused
	 * @param s0              starting path length in meters
	 * @param sMax            maximal path length in meters
	 * @param rho0            starting rho
	 */
	public void set(double s0, double sMax, double rho0) {
		_s0 = s0;		
		_totS = 0;
		_prevS = 0;

		_sMax = sMax;
		_startSign = sign(rho0);
	}
	
	@Override
	public boolean stopIntegration(double s, double[] u) {
		
		double currentRho = Math.hypot(u[0], u[1]);
		_totS = s;

		// within accuracy?
		//note this could also result with s > smax
		if (Math.abs(currentRho - _targetRho) < _accuracy) {
            copy(u, _uf);
			return true;
		}

		//stop (and backup/reset to prev) if we crossed the boundary or exceeded smax
		if ((getFinalS() > _sMax) || (sign(currentRho) != _startSign)) {
			_totS = _prevS;
			return true;
		}
				
		//copy current to previous
		_prevS = _totS;
        copy(u, _uf);
		return false;
	}

	/**
	 * Get the final path length in meters
	 * 
	 * @return the final path length in meters
	 */
	@Override
	public double getFinalS() {
		return _s0 + _totS;
	}
	
	/**
	 * Get the final value of the state vector
	 * @return
	 */
	public double[] getFinalU() {
		if (_uf == null) {
			System.err.println("Returning null final u");
		}
		return _uf;
	}

	//get the sign based on the current rho
	private int sign(double currentRho) {
		return ((currentRho < _targetRho) ? -1 : 1);
	}
	
	//array copy for state vectors
	private void copy(double src[], double[] dest) {
		System.arraycopy(src, 0, dest, 0, _dim);
	}



}