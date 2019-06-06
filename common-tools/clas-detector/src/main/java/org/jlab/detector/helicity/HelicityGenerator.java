package org.jlab.detector.helicity;

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Helicity Pseudo-Random Sequence.
 *
 * This calculates the first helicity state in each pattern (quartet/octet/toggle),
 * according to the specs here:
 * 
 * https://hallaweb.jlab.org/equipment/daq/HelicityUsersGuideFeb4.pdf
 *
 * Note, the convention above, and internally in this class, is HIGH=0=+ and LOW=1=-.
 *
 * Note, this pseudorandom generator operates on raw helicity, and does not produce
 * the same sequence if you just reverse all the states, so the HWP-corrected
 * helicity cannot be used here.
 * 
 * User calls addState() until initialized()==true, then getState().
 *
 * Public interface requires HelicityBit (not integers).
 *
 * @author baltzell
 */
public final class HelicityGenerator {

    public static final int REGISTER_SIZE=30;
    private final List<Integer> states=new ArrayList<>();
    private int register=0;
    private long timestamp=0;

    public HelicityGenerator(){}

    /**
     * Get the next bit in the sequence.
     * Requires initialized()==true.
     */
    private int getNextBit() {
        final int bit7  = (this.register>>6)  & 0x1;
        final int bit28 = (this.register>>27) & 0x1;
        final int bit29 = (this.register>>28) & 0x1;
        final int bit30 = (this.register>>29) & 0x1;
        final int nextBit = (bit30^bit29^bit28^bit7) & 0x1;
        return nextBit;
    }

    /**
     * Get the timestamp of the first state in the generator sequence.
     * @return timestamp (4ns)
     */
    public long getTimestamp() {
        return timestamp;
    }
    
    /**
     * Shift the register with the next state.
     * Requires initialized()==false.
     */
    private void shiftRegister(int state) {
        this.states.add(state);
        this.register = ( state |  (this.register<<1) ) & 0x3FFFFFFF;
    }

    /**
     * Shift the register with the next state.
     * Requires initialized()==true.
     */
    private void shiftRegister() {
        this.shiftRegister(this.getNextBit());
    }

    /**
     * Reset the generator, clearing all added states, in order to
     * reuse this {@link HelicityGenerator} for a new sequence.
     */
    public void reset() {
        this.states.clear();
        this.register=0;
        this.timestamp=-1;
    }

    /**
     * Test whether the generator is sufficiently initialized such that
     * {@link getState(int)} method can be called, based on whether the
     * number of added states is at least {@link REGISTER_SIZE}.
     *
     * @return whether the sequence is initialized
     */
    public boolean initialized() {
        return this.states.size() >= REGISTER_SIZE;
    }

    /**
     * Get the number of states currently in the generator's sequence.
     * 
     * @return size 
     */
    public int size() {
        return this.states.size();
    }
    
    /**
     * Let the user add a state to initialize the sequence.
     * 
     * This must be the first helicity state in the next pattern.  States
     * must be added serially and without skipped patterns.  Requires 
     * {@link initialized} is false and a defined HelicityBit, else throws
     * an exception..
     *
     * @param bit the HelicityBit to add to the sequence.  This must
     * be the raw helicity, e.g. HelicityState.getHelicityRaw(), not the
     * HWP-corrected version.
     */
    public void addState(HelicityBit bit) {
        if (bit==HelicityBit.UDF)
            throw new RuntimeException("Helicity state undefined.");
        if (this.initialized())
            throw new RuntimeException("Already initialized");
        this.shiftRegister(bit==HelicityBit.PLUS?0:1);
    }

    /**
     * Let the user add a state to initialize the sequence.
     * 
     * This just calls {@link #addState(HelicityBit)} with the raw
     * {@link #HelicityBit} from the given {@link #HelicityState}.
     *
     * This must be the first helicity state in the next pattern.  States
     * must be added serially and without skipped patterns.  Requires
     * {@link initialized} is false and a defined {@link #HelicityState}.
     * 
     * @param state the HelicityState to add to the sequence.
     */
    public void addState(HelicityState state) {
        this.addState(state.getHelicityRaw());
    }
   
    /**
     * Get the nth state in the sequence.
     * Requires initialized()==true.
     *
     * @param n number of states after the first one.
     * @return the nth HelicityBit in the sequence.
     */
    public HelicityBit getState(final int n) {
        if (!this.initialized())
            throw new RuntimeException("Not initialized.");
        if (n < 0)
            throw new RuntimeException("Invalid argument.");
        while (this.states.size() <= n)
            this.shiftRegister();
        return this.states.get(n) == 0 ?
            HelicityBit.PLUS : HelicityBit.MINUS;
    }
    
   
    /**
     * Initialize with a list of states.
     * 
     * The states are first time-ordered, and error-checking is done to find the
     * first valid sequence of suffecient length in the list and use it to
     * initialize the generator, otherwise the return value will be false.
     * @param states list of HelicityState objects
     * @return success of initializing the generator 
     */
    public final boolean initialize(List<HelicityState> states) {

        // make sure they're time-ordered:
        Collections.sort(states);
        
        // reset this generator:
        this.reset();
      
        // these will be the first valid sequence found in the input states:
        List<HelicityState> genStates=new ArrayList<>();

        for (int ii=0; ii<this.states.size(); ii++) {
            
            // any initial state will do:
            if (genStates.isEmpty()) {
                genStates.add(states.get(ii));
            }
            
            // bad pair sync, reset the sequence:
            else if (states.get(ii).getPairSync().value() ==
                genStates.get(genStates.size()-1).getPairSync().value()) {
                genStates.clear();
            }

            // bad pattern sync, reset the sequence:
            else if (genStates.size() > 2) {
                if (states.get(ii).getPatternSync().value()+
                    genStates.get(genStates.size()-1).getPatternSync().value()+
                    genStates.get(genStates.size()-2).getPatternSync().value()+
                    genStates.get(genStates.size()-3).getPatternSync().value() != 2) {
                    genStates.clear();
                }
            }
           
            else {

                // get time difference between states:
                final double seconds = (states.get(ii).getTimestamp() -
                        genStates.get(genStates.size()-1).getTimestamp()) / HelicitySequence.TIMESTAMP_CLOCK;
                
                // bad timestamp delta, reset the sequence:
                if (seconds < (1.0-0.5)/HelicitySequence.HELICITY_CLOCK ||
                    seconds > (1.0+0.5)/HelicitySequence.HELICITY_CLOCK) {
                    genStates.clear();
                }

                // passed all checks, add the state:
                else {
                    genStates.add(states.get(ii));
                }
            }
       
            // ok, we got enough to initialize the generator:
            // FIXME: we assume quartet here (the "4")
            if (genStates.size() >= REGISTER_SIZE*4+1) {
                
                int idxOffset = -1;

                // store all states timestamps to apply a modulo correction
                List <Double> timestamps = new ArrayList<>();
               
                for (int jj=0; jj<genStates.size(); jj++) {

                    // we've got enough valid, consecutive states, so correct
                    // the first state's timestamp and stop the investigation:
                    if (this.initialized()) {
                        this.timestamp=0;
                        for (int kk=0; kk<timestamps.size(); kk++) {
                            this.timestamp += timestamps.get(kk);
                        }
                        this.timestamp /= timestamps.size();
                        break;
                    }

                    if (genStates.get(jj).getPatternSync() == HelicityBit.MINUS) {
                        if (this.size() == 0) {
                            idxOffset = jj;
                        }
                        this.addState(states.get(ii));
                    }

                    // subtract off the nominal flip period:
                    if (this.size() > 0) {
                        timestamps.add(genStates.get(jj).getTimestamp() - (jj-idxOffset)/HelicitySequence.HELICITY_CLOCK*HelicitySequence.TIMESTAMP_CLOCK);
                    }
                }
            }
        }

        if (!this.initialized()) {
            this.reset();
        }

        return this.initialized();
    }
}
