package cnuphys.bCNU.graphics.toolbar;

import java.awt.event.MouseEvent;

import cnuphys.bCNU.component.MagnifyWindow;
import cnuphys.bCNU.graphics.container.IContainer;
import cnuphys.bCNU.view.BaseView;

/**
 * Used to rubber band a zoom.
 * 
 * @author heddle
 * 
 */
@SuppressWarnings("serial")
public class MagnifyButton extends ToolBarToggleButton {


    /**
     * Create the button for magnification
     * 
     * @param container
     *            the owner container.
     */
    public MagnifyButton(IContainer container) {
	super(container, "images/magnify.png", "Magnification");
	customCursorImageFile = "images/box_zoomcursor.gif";
    }
    
    /**
     * Handle a mouse enter (into the container) event (if this tool is active).
     * 
     * @param mouseEvent
     *            the causal event.
     */
    @Override
    public void mouseEntered(MouseEvent mouseEvent) {
	System.err.println("Mouse enter with mag glass");
    }
   
    
    /**
     * Handle a mouse exit (into the container) event (if this tool is active).
     * 
     * @param mouseEvent
     *            the causal event.
     */
    @Override
    public void mouseExited(MouseEvent mouseEvent) {
	System.err.println("Mouse exited with mag glass");
	MagnifyWindow.closeMagnifyWindow();
    }


    /**
     * Handle a mouse press (into the container) event (if this tool is active).
     * 
     * @param mouseEvent
     *            the causal event.
     */
    @Override
    public void mousePressed(MouseEvent mouseEvent) {
	System.err.println("Mouse pressed with mag glass");
    }
    

    /**
     * Handle a mouse move (into the container) event (if this tool is active).
     * 
     * @param mouseEvent
     *            the causal event.
     */
    @Override
    public void mouseMoved(MouseEvent mouseEvent) {
	System.err.println("Mouse moved with mag glass");
	BaseView view = container.getView();
	
	if (view == null) {
	    return;
	}
	
	view.handleMagnify(mouseEvent);
    }
   
    
}