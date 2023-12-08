%macro deleteAllMacroVars();
	%local macro_vars_to_delete;
 
	proc sql noprint;
		select name into :macro_vars_to_delete separated by " " 
		from dictionary.macros
		where scope="GLOBAL" 
		and not substr(name,1,3) =  "SYS"
		and not name = "GRAPHTERM";
	quit;
 
	%symdel &macro_vars_to_delete.;
%mend deleteAllMacroVars;
