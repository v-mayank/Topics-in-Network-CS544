# Topics-in-Network-CS544
Assignment to model various queueing processes(INQ,KOUTQ,iSLIP) in switches.
194101031.zip folder contains:-
1.Source code (ASSG2_19410103.cpp)
2.Report file (Report_ASSG2_194101031.pdf)
3.README.txt.

How to run source code:-
	1. Compile using g++ command.
	2. Run with .out file.(By default INQ will run).
	3. To run ISLIP pass "Q_ISLIP" argument, similerly "Q_KOUQ" for knockout and "Q_INQ" for INQ.
	4. You have to pass letter then 'underscore' then value.	
	5. To change default values pass arguments like:
		N_8   	for number of ports, 
		P_0.6 	for generation probability,
		K_0.6	for knockout value,
		Q_INQ 	for queue type,
		B_4	for buffer size &
		T_1000	for maximum time slot.
	
	6. Any order of arguments will run the code.
