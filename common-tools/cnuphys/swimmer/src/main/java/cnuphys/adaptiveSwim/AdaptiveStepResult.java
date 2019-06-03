package cnuphys.adaptiveSwim;

public class AdaptiveStepResult {
	
	//the stepsize we actually used
	private double _hUsed;
	
	//new stepsize
	private double _hNew;
	
	//new value of independent variable
	private double _sNew;

	public double getHNew() {
		return _hNew;
	}

	public void setHNew(double hNew) {
		_hNew = hNew;
	}

	public double getSNew() {
		return _sNew;
	}

	public void setSNew(double sNew) {
		_sNew = sNew;
	}
	
	public double getHUsed() {
		return _hUsed;
	}

	public void setHUsed(double hUsed) {
		_hUsed = hUsed;
	}


}
