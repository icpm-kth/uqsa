TITLE AKAP79
COMMENT
	automatically generated from an SBtab file
	date: Fri Oct 20 16:40:07 2023
ENDCOMMENT
NEURON {
	SUFFIX AKAP79 : OR perhaps POINT_PROCESS ?
	RANGE b_AKAP : input
	RANGE AKAR4pOUT : output
	RANGE kf_RiiP_cAMP_CaN__CaNXRii_cAMP : assigned
	RANGE kb_RiiPxCaN__RiiP_CaN : assigned
	RANGE kf_RiiP_CaN__RiixCaN : assigned
	RANGE kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN : assigned
	RANGE kf_RiiPxCaN__RiiP_CaN : assigned
	RANGE kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN : assigned
	RANGE kb_RiiP_CxcAMP__RiiP_C_cAMP : assigned
	RANGE Rii : compound
	RANGE cAMP : compound
	RANGE RiiP : compound
	RANGE Rii_C : compound
	RANGE RiiP_cAMP : compound
	RANGE RiiP_C : compound
	RANGE RiiP_C_cAMP : compound
	RANGE C : compound
	RANGE Rii_cAMP : compound
	RANGE Rii_C_cAMP : compound
	RANGE CaN : compound
	RANGE RiiP_CaN : compound
	RANGE RiiP_cAMP_CaN : compound
	RANGE AKAR4 : compound
	RANGE AKAR4_C : compound
	RANGE AKAR4p : compound
: USEION ca READ cai VALENCE 2 : sth. like this may be needed for ions you have in your model
}
CONSTANT {
}
PARAMETER {
	kf_Rii_C__RiiP_C = 0.033 (/millisecond): a kinetic parameter
	kf_RiiP_CxcAMP__RiiP_C_cAMP = 0.000496 (/micromolarity-millisecond): a kinetic parameter
	kf_RiiP_cAMPxC__RiiP_C_cAMP = 5.45e-06 (/micromolarity-millisecond): a kinetic parameter
	kb_RiiP_cAMPxC__RiiP_C_cAMP = 1.56e-05 (/millisecond): a kinetic parameter
	kb_RiiPXcAMP__RiiP_cAMP = 1.6e-06 (/millisecond): a kinetic parameter
	kf_RiiPXcAMP__RiiP_cAMP = 1.5e-05 (/micromolarity-millisecond): a kinetic parameter
	kf_RiiPxC__RiiP_C = 3.8e-05 (/micromolarity-millisecond): a kinetic parameter
	kb_RiiPxC__RiiP_C = 2.6e-06 (/millisecond): a kinetic parameter
	kf_cAMPxRii__Rii_cAMP = 1.5e-05 (/micromolarity-millisecond): a kinetic parameter
	kb_cAMPxRii__Rii_cAMP = 1.6e-06 (/millisecond): a kinetic parameter
	kf_Rii_CxcAMP__Rii_C_cAMP = 0.000496 (/micromolarity-millisecond): a kinetic parameter
	kb_Rii_CxcAMP__Rii_C_cAMP = 0.001413 (/millisecond): a kinetic parameter
	kf_RiixC__Rii_C = 0.0021 (/micromolarity-millisecond): a kinetic parameter
	kf_Rii_cAMPxC__Rii_C_cAMP = 0.0002984 (/micromolarity-millisecond): a kinetic parameter
	kb_Rii_cAMPxC__Rii_C_cAMP = 1.8e-05 (/millisecond): a kinetic parameter
	kf_Rii_C_cAMP__RiiP_C_cAMP = 0.033 (/millisecond): a kinetic parameter
	kb_RiixC__Rii_C = 3e-07 (/millisecond): a kinetic parameter
	AKAPoff_1 = 0.0026 (/millisecond): a kinetic parameter
	AKAPoff_3 = 0.02 (/millisecond): a kinetic parameter
	AKAPon_1 = 0.00045 (/millisecond): a kinetic parameter
	AKAPon_3 = 0.002 (/millisecond): a kinetic parameter
	kf_C_AKAR4 = 1.8e-05 (/micromolarity-millisecond): a kinetic parameter
	kb_C_AKAR4 = 0.000106 (/millisecond): a kinetic parameter
	kcat_AKARp = 0.0102 (/millisecond): a kinetic parameter
	kmOFF = 100 (micromolarity): a kinetic parameter
	kmON = 1 (micromolarity): a kinetic parameter
	KD_T = 0.7 (micromolarity): a kinetic parameter
	b_AKAP  = 0 (1) : an input
	AKAR4_ConservedConst = 0.2 : the total amount of a conserved sub-set of states
	CaN_ConservedConst = 1.5 : the total amount of a conserved sub-set of states
	Rii_C_ConservedConst = 0.63 : the total amount of a conserved sub-set of states
	cAMP_ConservedConst = 0 : the total amount of a conserved sub-set of states
	Rii_ConservedConst = 6.3 : the total amount of a conserved sub-set of states
}
ASSIGNED {
	time (millisecond) : alias for t
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP : a pre-defined algebraic expression
	kb_RiiPxCaN__RiiP_CaN : a pre-defined algebraic expression
	kf_RiiP_CaN__RiixCaN : a pre-defined algebraic expression
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN : a pre-defined algebraic expression
	kf_RiiPxCaN__RiiP_CaN : a pre-defined algebraic expression
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN : a pre-defined algebraic expression
	kb_RiiP_CxcAMP__RiiP_C_cAMP : a pre-defined algebraic expression
	reaction_51 : a flux, for use in DERIVATIVE mechanism
	reaction_14 : a flux, for use in DERIVATIVE mechanism
	reaction_12 : a flux, for use in DERIVATIVE mechanism
	reaction_43 : a flux, for use in DERIVATIVE mechanism
	reaction_23 : a flux, for use in DERIVATIVE mechanism
	reaction_78 : a flux, for use in DERIVATIVE mechanism
	reaction_56 : a flux, for use in DERIVATIVE mechanism
	reaction_76 : a flux, for use in DERIVATIVE mechanism
	reaction_62 : a flux, for use in DERIVATIVE mechanism
	reaction_58 : a flux, for use in DERIVATIVE mechanism
	reaction_44 : a flux, for use in DERIVATIVE mechanism
	reaction_33 : a flux, for use in DERIVATIVE mechanism
	reaction_48 : a flux, for use in DERIVATIVE mechanism
	reaction_37 : a flux, for use in DERIVATIVE mechanism
	reaction_1 : a flux, for use in DERIVATIVE mechanism
	reaction_2 : a flux, for use in DERIVATIVE mechanism
	AKAR4 : computed from conservation law
	CaN : computed from conservation law
	Rii_C : computed from conservation law
	cAMP : computed from conservation law
	Rii : computed from conservation law
	AKAR4pOUT : an observable
}
PROCEDURE assign_calculated_values() {
	time = t : an alias for the time variable, if needed.
	AKAR4 = AKAR4_ConservedConst - (AKAR4_C+AKAR4p) : conservation law
	CaN = CaN_ConservedConst - (RiiP_CaN+RiiP_cAMP_CaN) : conservation law
	Rii_C = Rii_C_ConservedConst - (RiiP_C+RiiP_C_cAMP+C+Rii_C_cAMP+AKAR4_C) : conservation law
	cAMP = cAMP_ConservedConst - (RiiP_cAMP+RiiP_C_cAMP+Rii_cAMP+Rii_C_cAMP+RiiP_cAMP_CaN) : conservation law
	Rii = Rii_ConservedConst - (RiiP+RiiP_cAMP+Rii_cAMP+RiiP_CaN+RiiP_cAMP_CaN-C-AKAR4_C) : conservation law
	kf_RiiP_cAMP_CaN__CaNXRii_cAMP = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1 : assignment for expression kf_RiiP_cAMP_CaN__CaNXRii_cAMP
	kb_RiiPxCaN__RiiP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3 : assignment for expression kb_RiiPxCaN__RiiP_CaN
	kf_RiiP_CaN__RiixCaN = b_AKAP * AKAPon_1 + (1 - b_AKAP) * AKAPoff_1 : assignment for expression kf_RiiP_CaN__RiixCaN
	kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP*AKAPon_3  +  (1 - b_AKAP)* AKAPoff_3 : assignment for expression kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN
	kf_RiiPxCaN__RiiP_CaN = b_AKAP * ((kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP)/kmON ) + (1 - b_AKAP) * (kb_RiiPxCaN__RiiP_CaN + kf_RiiP_cAMP_CaN__CaNXRii_cAMP) / kmOFF : assignment for expression kf_RiiPxCaN__RiiP_CaN
	kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN = b_AKAP * ((kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmON) + (1 - b_AKAP) * (kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN + kf_RiiP_CaN__RiixCaN)/kmOFF : assignment for expression kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN
	kb_RiiP_CxcAMP__RiiP_C_cAMP = kf_RiiP_CxcAMP__RiiP_C_cAMP * KD_T : assignment for expression kb_RiiP_CxcAMP__RiiP_C_cAMP
	reaction_51 = kf_Rii_C__RiiP_C*Rii_C : flux expression reaction_51
	reaction_14 = kf_RiiPxC__RiiP_C*RiiP*C - kb_RiiPxC__RiiP_C*RiiP_C : flux expression reaction_14
	reaction_12 = kf_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C*cAMP - kb_RiiP_CxcAMP__RiiP_C_cAMP*RiiP_C_cAMP : flux expression reaction_12
	reaction_43 = kf_RiiPXcAMP__RiiP_cAMP*cAMP*RiiP - kb_RiiPXcAMP__RiiP_cAMP*RiiP_cAMP : flux expression reaction_43
	reaction_23 = kf_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_cAMP*C - kb_RiiP_cAMPxC__RiiP_C_cAMP*RiiP_C_cAMP : flux expression reaction_23
	reaction_78 = kf_cAMPxRii__Rii_cAMP*cAMP*Rii - kb_cAMPxRii__Rii_cAMP*Rii_cAMP : flux expression reaction_78
	reaction_56 = kf_Rii_CxcAMP__Rii_C_cAMP*Rii_C*cAMP - kb_Rii_CxcAMP__Rii_C_cAMP*Rii_C_cAMP : flux expression reaction_56
	reaction_76 = kf_Rii_cAMPxC__Rii_C_cAMP*Rii_cAMP*C - kb_Rii_cAMPxC__Rii_C_cAMP*Rii_C_cAMP : flux expression reaction_76
	reaction_62 = kf_Rii_C_cAMP__RiiP_C_cAMP*Rii_C_cAMP : flux expression reaction_62
	reaction_58 = kf_RiixC__Rii_C*Rii*C - kb_RiixC__Rii_C*Rii_C : flux expression reaction_58
	reaction_44 = kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN : flux expression reaction_44
	reaction_33 = kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN : flux expression reaction_33
	reaction_48 = kf_RiiP_CaN__RiixCaN*RiiP_CaN : flux expression reaction_48
	reaction_37 = kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN : flux expression reaction_37
	reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C : flux expression reaction_1
	reaction_2 = kcat_AKARp*AKAR4_C : flux expression reaction_2
}
STATE {
	: Rii is calculated via Conservation Law
	: cAMP is calculated via Conservation Law
	RiiP (micromole/liter) : a state variable
	: Rii_C is calculated via Conservation Law
	RiiP_cAMP (micromole/liter) : a state variable
	RiiP_C (micromole/liter) : a state variable
	RiiP_C_cAMP (micromole/liter) : a state variable
	C (micromole/liter) : a state variable
	Rii_cAMP (micromole/liter) : a state variable
	Rii_C_cAMP (micromole/liter) : a state variable
	: CaN is calculated via Conservation Law
	RiiP_CaN (micromole/liter) : a state variable
	RiiP_cAMP_CaN (micromole/liter) : a state variable
	: AKAR4 is calculated via Conservation Law
	AKAR4_C (micromole/liter) : a state variable
	AKAR4p (micromole/liter) : a state variable
}
INITIAL {
	: Rii cannot have initial values as it is determined by conservation law
	: cAMP cannot have initial values as it is determined by conservation law
	 RiiP = 0 : initial condition
	: Rii_C cannot have initial values as it is determined by conservation law
	 RiiP_cAMP = 0 : initial condition
	 RiiP_C = 0 : initial condition
	 RiiP_C_cAMP = 0 : initial condition
	 C = 0 : initial condition
	 Rii_cAMP = 0 : initial condition
	 Rii_C_cAMP = 0 : initial condition
	: CaN cannot have initial values as it is determined by conservation law
	 RiiP_CaN = 0 : initial condition
	 RiiP_cAMP_CaN = 0 : initial condition
	: AKAR4 cannot have initial values as it is determined by conservation law
	 AKAR4_C = 0 : initial condition
	 AKAR4p = 0 : initial condition
}
BREAKPOINT {
	SOLVE ode METHOD cnexp
	assign_calculated_values() : procedure
}
DERIVATIVE ode {
	: Compound Rii with initial condition 6.3 had derivative -reaction_78-reaction_58+reaction_48, but is calculated by conservation law.
	: Compound cAMP with initial condition 0 had derivative -reaction_12-reaction_43-reaction_78-reaction_56, but is calculated by conservation law.
	RiiP' = -reaction_14-reaction_43-reaction_44 : affects compound RiiP
	: Compound Rii_C with initial condition 0.63 had derivative -reaction_51-reaction_56+reaction_58, but is calculated by conservation law.
	RiiP_cAMP' = reaction_43-reaction_23-reaction_33 : affects compound RiiP_cAMP
	RiiP_C' = reaction_51+reaction_14-reaction_12 : affects compound RiiP_C
	RiiP_C_cAMP' = reaction_12+reaction_23+reaction_62 : affects compound RiiP_C_cAMP
	C' = -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2 : affects compound C
	Rii_cAMP' = reaction_78-reaction_76+reaction_37 : affects compound Rii_cAMP
	Rii_C_cAMP' = reaction_56+reaction_76-reaction_62 : affects compound Rii_C_cAMP
	: Compound CaN with initial condition 1.5 had derivative -reaction_44-reaction_33+reaction_48+reaction_37, but is calculated by conservation law.
	RiiP_CaN' = reaction_44-reaction_48 : affects compound RiiP_CaN
	RiiP_cAMP_CaN' = reaction_33-reaction_37 : affects compound RiiP_cAMP_CaN
	: Compound AKAR4 with initial condition 0.2 had derivative -reaction_1, but is calculated by conservation law.
	AKAR4_C' = reaction_1-reaction_2 : affects compound AKAR4_C
	AKAR4p' = reaction_2 : affects compound AKAR4p
}
PROCEDURE observables_func() {
	AKAR4pOUT = (AKAR4p*5)*71.67+100 : Output ID AKAR4pOUT
}
