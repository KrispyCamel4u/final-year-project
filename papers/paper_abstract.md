 
# Low-Latency distance protective relay on fpga

	-the fpga offer the high parallelism and throughput with minimum latency
### option 1 ###
	- the very first block is the DFT(discrete fourier transform) 
	- DFT is the sampled version of DTFT(discreate time fourier transform)
	- there is a CORDIC unit which provides the trignometric functions
	- then the overcurrent relay logic that is the fault detection logic
	- then protection element and trip logic wiich deteremines the impedance and based on that it takes the decision finally it issues the trip signal
### option 2 ###
	-- in parallel the instantaneous-signal based distance protection block which can also issues the trip signal 	

### each module specs ###

#### DFT ####
	- DFT has improved harmonic immunity
	- need dc offset removel
	- for DFT one complete cycle and two samples for dc offset removel 
	- well there is some modelling of dc offset with some decaying functions which i don't understand
	- a FSM is there for which detects from the previous cycle weather the fault is detected or not if detected than it will send the offset signal to include the offset part in the DFT which a liitle unclear to me

#### CORDIC BLOCK
	-- for providing the non-linear and trignometric functions
	-- better precison than lookup tables method
	-- rotation and vector mode two ways to impliment
	-- piplined algo on fpga
	-- (<CORDIC algo based on fixed point data>)

#### fault detection
	-- based on overcurrent method
	-- FCDFT have drawbacks such as sensitivity to even harmonics and decaying dc signal
	-- using the strategy such as when the current exceeds the threshold three times consecutively then it is treated as the abnormal condition and fault detection is sent out

#### distance protection element
	-- three phase to  ground and three phase to phase relays
	-- (<operating vector and polarising vector>)
	-- if the angle between Vo and Vp is less than or equal to 90 than it is in the zone and otherwise outside the zone 
	-- different zones can be achieved by varring the diameter of the circle

#### instantneous signal based distance relay element
	-- high frequency signals are filtered out
	---
## fpga
	-- CLB -configurable logic block
	-- IOB : input output block
	-- PIB : programable interconnection block
	-- PSCAD/EMTDC (some software used to simulate faults)

## test setup
	-- two syncrhonouys generator
	-- 210KM long transmission
	-- latencies 2.02us and 0.35us with clock 100MHz
	-- DFT module consumes largest latency
	-- instantaneous took 35 cyles

## hardware experimental results
	-- testing under 4 types of fault condition
	-- whenever the faults occurs the transients are invoked 
	-- before fault the impedance will be outside the mho criteria circle but after the fault it moves within the operating circle
	-- increasing sampling time can shorten the time of detection
	
#### fpga advantages
		-- parallel architecture
		-- each task has its own alocated hardware
		-- also paralled distributed memory can be given in fpga
		-- however in DSP based systems the algorithm relies on the to from information from the memory and also the the resourses are shared hence a small mistake can break down whole system
		-- hence fpga offers higher reliability
		-- also troubleshooting is easy becouse of independent seperate functional blocks can be tested
		-- since the operating clock speed and the hardware resource usages are low as compared to the DSPs which may operated in around 1GHz which inturn causes the heat and life cylce issue on the other hand fpga doen't have this problem
		-- partial reconfiguration (PR) features in which partial configuration bit file can modify the operating fpga which gives the ability of on-site programming.
	

### conclusion
	-- with all the cases the it operated fairly well with low latencies and high pipelined nature of fpga. 
	-- can be used to inhance the the functionality of the module to make it operate as a smart meter or anything else







# paper 2   

For transformer differential protection is employed in general
	-- but the CTs suffers from the saturation problem
	-- subjected to false trapping due to the inrush current :since the inrush current is rich in second order harmonics(not able to relate much)

	-- the second harmonins in inrush current in transmorfersm is very less which are connected to the long transmission lines

	-- there was a factor to account for the second harmonics but since in case of the long transmission lines the second harmonics is less which can resul t in muloperation in the relay which depend on the second harmonic restraint
	-- the persence of nonlinear devices and the underground cables(<doubt>) the fault current has increased ammount of second haronics and becouse of the restrain the tripping may get delayed.
	
### the present work
	-- have two parts:
		disturbance detection and fault discrimination
	(<what is the first level high-frequency>)
	(<wavelet tranform>)
	-- directional signal is obtained by the product frequency details of current and voltage
	(<EMTP/ATP>)
	(<bolted faults>)
	(<winding model didn't understand>)
	
	-- transformer faults have low frequency oscillatory transients <5Khz
	-- wavelet transform is used for the extraticg the information from transients
	-- WT(waveler transform) can give time and frequency info simultaneously
	-- fast fourier transform(any efficient algo calculating DFT) works well for smooth and uniform frequencies 
	-- WT works better with signals with  sharp discontinuities and non stationary nature
	-- WT less commutations (O(N)) ,FFT (O(NlogN)) N being the number of samples
	-- mother wavelet is a functions used for detecting and localising(<>) different types of faults transients,"daubechies wavelet fucntion works for the protection ckt
	-- (<anti-aliasing filters>) appers to be same as the LPF
	-- (<down sampled>)
	-- (<high frequency detail coefficents>) in the range of 500-1000Hzare obtained 
	-- (<first difference of a sine wave>)


### signal preprocessing
	-- for 50Hz system frequency band of 0-1000Hz is more informative 
	-- a sampling frequency of 2KHz can be chosen
	-- 1. transducesrs and isolations 2. anti aliasing filter 3. sample and hold 4. ADC 5. multiplexer
	-- a low pass filter is used to band limit the signal to the half of sampling frequency
	-- a butterworth filter is good as to maintain the balance between the time domain response and the frequency domain response 
	-- (<as the frequency responce is sharper the time response becomes worse>)

### disturbance detection
	
	-- the disturbance is detected with the phase voltages as the change in voltage is instantaneous rather than current
	-- the detect signal goes high if distrubance is detected in any of the signal
	-- the thresold depends on expected voltage level, frequency and wavelet

### fault discrimination
	
	-- internal and external signals are decided based on the power directionl signal
	-- direction power signals are the cumulative sum of the the phases of for each side (this is done for high frequecies) (multiplication of high frequency detail signal of current and voltages give the hfd power signal) it is negligible when no foault occurs as the high frequency signal is absent
	-- the fault intrduces the high frequency oscillations and the high frequency power flows
	-- after 5ms the direction of Phv and Plv is determined ,if both have same direction then internal fault otherwise externel
	-- RTWT real-time window target of MATLAB	
	-- how will it work in case of ct saturations
	



# paper 3 high speed digital distance relaying scheme
	
	-- FCDFT requires the full cycle samples
	-- for faster speed HCDFT(half cycle DFT) but losed the ability to reject the even harmonics, also dificult to remove the DC offset
	-- usign Wavelet transformation multiresolution analysis fundamental phasor can be ectracted faster than the FCDFT , but performance is affected by the ddc
	-- in this paper phaselet-based distance scheme is developed

### phaselet based distance protective relay
	 
	-- the FCDFT has a fixed window size that spans between two partial window during the fault conditions to address this  issur phaselet based variable window filtering techique is developed
	-- phaselets are the partial integrals of the product of the waveform samples and the sinosoidal and cosinosoidal coefficients(>)
	-- under normal conditions the window size fixed to N when fault occurs the window will shrink to the size of P(number of samples per phaselets)
	-- the FCDFT have zero gains for the other harmonics and as in the phaselet filtering window increases the harmonic rejection increases and dc offset also
	-- a disurbance ditector is used for the purpose of initiating the adaptive phasor estimation for this previous two cycle are taken into account
	-- the change in magnitude is evaluated against 2 times the cutoff value

#### transient monitoring and error analysis
	
	-- the estimation error are large when the window size is less the the 40% of the full cycle which makes the initial phasor estimation unusable
	-- phase angle errors are very important and has significant effect
	-- a trip delay routine is used to prevent false trips due to phase angle errors

#### adaptive mho characteristics
	
	-- when the window size is less than 40% the reach is set to 0
	-- close in fault which occurs near to the relay
	-- the close in faults can be detected quickly with less window size 
	
#### trip delay routine
	
	-- activated when the impedace trajectory enter the mho circle for first time
	-- for the couple or more phaselet are recorded and then the trip decision is made otherwise the routine is reset and wait for the new activation
	-- relay accuracy is affected by the voltage level

### hardware implimentation

	-- mimic filter , adaptive phasor estimation, 3 phase to symmetrical components, impedance calculation, adaptive mho characteistics
	
### hardware in  the looop testing and results

	-- a 12 bus ieee standard test system 

#### operating time
	
	-- noted for different system for different inception angle 
	-- for line ground faults the phaselet method performed better than the other methods ,for three phase faults the operating speed is faster , for line to line is slower but within .85 cycles    but faster than the FCDFT 
	-- SIR(source impedance ratio)

#### reliability
	
	-- with higher SIR a longer trip delay is needed
	-- 2% possibility of false trip in ground faults
	-- 1% for the line to line

#### effect of CT saturation
	
	-- the ct satuaration doesn't occur instantly 
	-- in case of ct saturation also the for the half of the cycle it will be same as the non saturation case but than the impedance will increase and then settle to the same final value and the phaselet based relay takes the triping decision in the time when the ct satuatioin and without satuaration impedance signal overlaps so it the ct satuation doesn't effect it
	-- CCVT : coupling capacitor voltage transformer
	-- the ccvt transients errors results in smaller voltages and hence smaller apparent impedances
	-- the CCVT transients affects this relay and larger trip delay is needed

#### effect of energisation of shunt capacitor
	
	-- when the capacitor is energised or de-energiesed then a high freqquency high magnitude current flows through the capacitor
	-- the disturbance is picked up but the relay is no maloperated because the impedance doesn't come in the reach zone




# paper 4 embedding synchronized measurement technoloagy

	-- systainable development -> smartgrid: according to europian smart grid is as : electricity networks that can intelligently integrate the action of all users connected to it -generators, consumers and tohose that do both in  order to efficiently deliver sustainabke economic way of supply
	-- reliability : amount of time without interupption 
	-- power quality : clean voltage and current signals 
	-- DER : distributed energy resources 
	
### distribution automation system overview
		
		-- DAS enables electricity utilities to monitor cordinate  and operate distribution components in as real-time remote mode
		-- DMS distribution management system :software 















# paper 5 hormonic phasor estimator for p class
	
	-- high impedance fault occuring in distribution network causes harmonin signal in current signal
	-- Arcing faults can also introduce the harmonics in currents
	-- such harmonic info can be used in high impedance fault identification
	-- PMUs :P class for protection application and M class for monitoring and mesurement applications
	-- (in IEEE standards the reference algorithm for P class takes two cycles long)
	-- goal of this paper is two make two cylce harmoni PMU
	-- the faults currents may contains the decaying dc signal which make hard in achieving the the high accuracy and low complexity
	-- THE DFT is preferred for its simplicity for harmonic PMUs but is has large errors under frequency deviation, (<harmonic oscillations>), dc offsets
	-- DFT uses the static phasor model wheras the Taylor Fourier Transform(TFT) uses dynamic harmonic phasor based on Taylor signal model
	-- (<finite impulse response filter>) have higher passband and stopbands than the equivalent filters of DFT.
	-- the adaptive TFT can widen the passbands and stopbands by a (<first frequency estimation>)
	-- (<interharmonics>)
	-- for the interharmonics Taylor Fourier multifrequency model was proposed
	-- Fast TFM reduces the complixity of simplifying procedures of seeking frequency components
	--


