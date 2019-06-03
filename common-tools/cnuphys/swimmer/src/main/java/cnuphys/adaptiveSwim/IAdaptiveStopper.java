package cnuphys.adaptiveSwim;

public interface IAdaptiveStopper {

	/**
	 * Given the current state of the integration, should we stop? This allows the
	 * integration to stop, for example, if some distance from the origin has been
	 * exceeded or if the independent variable passes some threshold. It won't be
	 * precise, because the check may not happen on every step, but it should be
	 * close.
	 * 
	 * @param s the current value of the independent variable (typically pathlength)
	 * @param u the current state vector (typically [x, y, z, vx, vy, vz])
	 * @return <code>true</code> if we should stop now.
	 */
	public boolean stopIntegration(double s, double u[]);

	/**
	 * Get the final independent variable (typically path length in meters)
	 * 
	 * @return the final independent variable (typically path length in meters)
	 */
	public double getFinalS();


}
