TITLE AKAP79
COMMENT
	automatically generated from an SBtab file
	date: Mon Aug 15 16:34:20 2022
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
	kmOFF = 0.1 (/millisecond): a kinetic parameter
	kmON = 0.001 (/millisecond): a kinetic parameter
	KD_T = 0.0007 (/millisecond): a kinetic parameter
	b_AKAP  = 0 (1) : an input
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
	reaction_44_ : a flux, for use in DERIVATIVE mechanism
	reaction_33_ : a flux, for use in DERIVATIVE mechanism
	reaction_4_8 : a flux, for use in DERIVATIVE mechanism
	reaction_3_7 : a flux, for use in DERIVATIVE mechanism
	reaction_1 : a flux, for use in DERIVATIVE mechanism
	reaction_2 : a flux, for use in DERIVATIVE mechanism
	AKAR4pOUT : an observable
}
PROCEDURE assign_calculated_values() {
	time = t : an alias for the time variable, if needed.
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
	reaction_44_ = kf_RiiPxCaN__RiiP_CaN*RiiP*CaN - kb_RiiPxCaN__RiiP_CaN*RiiP_CaN : flux expression reaction_44
	reaction_33_ = kf_CaNxRiiP_cAMP__RiiP_cAMP_CaN*CaN*RiiP_cAMP - kb_CaNxRiiP_cAMP__RiiP_cAMP_CaN*RiiP_cAMP_CaN : flux expression reaction_33
	reaction_4_8 = kf_RiiP_CaN__RiixCaN*RiiP_CaN : flux expression reaction_48
	reaction_3_7 = kf_RiiP_cAMP_CaN__CaNXRii_cAMP*RiiP_cAMP_CaN : flux expression reaction_37
	reaction_1 = kf_C_AKAR4*C*AKAR4 - kb_C_AKAR4*AKAR4_C : flux expression reaction_1
	reaction_2 = kcat_AKARp*AKAR4_C : flux expression reaction_2
}
STATE {
	Rii (micromole/liter) : a state variable
	cAMP (micromole/liter) : a state variable
	RiiP (micromole/liter) : a state variable
	Rii_C (micromole/liter) : a state variable
	RiiP_cAMP (micromole/liter) : a state variable
	RiiP_C (micromole/liter) : a state variable
	RiiP_C_cAMP (micromole/liter) : a state variable
	C (micromole/liter) : a state variable
	Rii_cAMP (micromole/liter) : a state variable
	Rii_C_cAMP (micromole/liter) : a state variable
	CaN (micromole/liter) : a state variable
	RiiP_CaN (micromole/liter) : a state variable
	RiiP_cAMP_CaN (micromole/liter) : a state variable
	AKAR4 (micromole/liter) : a state variable
	AKAR4_C (micromole/liter) : a state variable
	AKAR4p (micromole/liter) : a state variable
}
INITIAL {
	 Rii = 6.3 : initial condition
	 cAMP = 0 : initial condition
	 RiiP = 0 : initial condition
	 Rii_C = 0.63 : initial condition
	 RiiP_cAMP = 0 : initial condition
	 RiiP_C = 0 : initial condition
	 RiiP_C_cAMP = 0 : initial condition
	 C = 0 : initial condition
	 Rii_cAMP = 0 : initial condition
	 Rii_C_cAMP = 0 : initial condition
	 CaN = 1.5 : initial condition
	 RiiP_CaN = 0 : initial condition
	 RiiP_cAMP_CaN = 0 : initial condition
	 AKAR4 = 0.2 : initial condition
	 AKAR4_C = 0 : initial condition
	 AKAR4p = 0 : initial condition
}
BREAKPOINT {
	SOLVE ode METHOD cnexp
	assign_calculated_values() : procedure
}
DERIVATIVE ode {
	Rii' = -reaction_78-reaction_58+reaction_4_8 : affects compound with ID Rii
	cAMP' = -reaction_12-reaction_43-reaction_78-reaction_56 : affects compound with ID cAMP
	RiiP' = -reaction_14-reaction_43-reaction_44_ : affects compound with ID RiiP
	Rii_C' = -reaction_51-reaction_56+reaction_58 : affects compound with ID Rii_C
	RiiP_cAMP' = reaction_43-reaction_23-reaction_33_ : affects compound with ID RiiP_cAMP
	RiiP_C' = reaction_51+reaction_14-reaction_12 : affects compound with ID RiiP_C
	RiiP_C_cAMP' = reaction_12+reaction_23+reaction_62 : affects compound with ID RiiP_C_cAMP
	C' = -reaction_14-reaction_23-reaction_76-reaction_58-reaction_1+reaction_2 : affects compound with ID C
	Rii_cAMP' = reaction_78-reaction_76+reaction_3_7 : affects compound with ID Rii_cAMP
	Rii_C_cAMP' = reaction_56+reaction_76-reaction_62 : affects compound with ID Rii_C_cAMP
	CaN' = -reaction_44_-reaction_33_+reaction_4_8+reaction_3_7 : affects compound with ID CaN
	RiiP_CaN' = reaction_44_-reaction_4_8 : affects compound with ID RiiP_CaN
	RiiP_cAMP_CaN' = reaction_33_-reaction_3_7 : affects compound with ID RiiP_cAMP_CaN
	AKAR4' = -reaction_1 : affects compound with ID AKAR4
	AKAR4_C' = reaction_1-reaction_2 : affects compound with ID AKAR4_C
	AKAR4p' = reaction_2 : affects compound with ID AKAR4p
}
PROCEDURE observables_func() {
	AKAR4pOUT = AKAR4p : Output ID AKAR4pOUT
}
