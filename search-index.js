var searchIndex = JSON.parse('{\
"pcon":{"doc":"Prompt COuNter, a short kmer counter.","t":[0,0,0,0,0,0,0,0,3,3,13,13,3,13,4,13,13,4,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,5,3,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,5,13,4,13,13,6,13,11,11,11,11,11,11,11,11,11,11,11,11,11,11,3,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,3,11,11,11,11,11,11,11,11,11,11,11,11,11,11,13,13,13,13,3,4,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11],"n":["cli","count","counter","dump","error","serialize","solid","spectrum","Command","Count","Count","Csv","Dump","Dump","DumpType","Pcon","Solid","SubCommand","abundance","abundance","augment_args","augment_args","augment_args","augment_args_for_update","augment_args_for_update","augment_args_for_update","augment_subcommands","augment_subcommands_for_update","borrow","borrow","borrow","borrow","borrow","borrow_mut","borrow_mut","borrow_mut","borrow_mut","borrow_mut","clone","clone_into","cmp","command","command_for_update","eq","fmt","fmt","fmt","fmt","fmt","from","from","from","from","from","from_arg_matches","from_arg_matches","from_arg_matches","from_arg_matches","from_arg_matches_mut","from_arg_matches_mut","from_arg_matches_mut","from_arg_matches_mut","group_id","group_id","group_id","has_subcommand","input","inputs","into","into","into","into","into","kmer_size","outputs","outputs","partial_cmp","quiet","record_buffer","subcommand","timestamp","to_owned","to_possible_value","try_from","try_from","try_from","try_from","try_from","try_into","try_into","try_into","try_into","try_into","type_id","type_id","type_id","type_id","type_id","update_from_arg_matches","update_from_arg_matches","update_from_arg_matches","update_from_arg_matches","update_from_arg_matches_mut","update_from_arg_matches_mut","update_from_arg_matches_mut","update_from_arg_matches_mut","value_variants","verbosity","count","Counter","borrow","borrow_mut","clone","clone_into","cmp","count_fasta","count_fasta","count_fasta","count_fasta","count_fasta","default","eq","fmt","from","from_stream","from_stream","from_stream","from_stream","from_stream","get","get","get","get","get","get_raw","hash","into","k","new","new","new","new","new","partial_cmp","raw","serialize","to_owned","try_from","try_into","type_id","dump","DumpTypeFromStr","Error","IO","Log","Result","TypeNotMatch","borrow","borrow_mut","fmt","fmt","from","from","from","into","provide","source","to_string","try_from","try_into","type_id","Serialize","borrow","borrow_mut","csv","csv","csv","csv","csv","from","into","new","pcon","pcon","pcon","pcon","pcon","solid","solid","solid","solid","solid","try_from","try_into","type_id","Solid","borrow","borrow_mut","from","from_count","get","get_canonic","into","k","new","set","set_canonic","try_from","try_into","type_id","FirstMinimum","PercentAtLeast","PercentAtMost","Rarefaction","Spectrum","ThresholdMethod","borrow","borrow","borrow_mut","borrow_mut","eq","fmt","from","from","from_count","get_threshold","into","into","try_from","try_from","try_into","try_into","type_id","type_id"],"q":["pcon","","","","","","","","pcon::cli","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","pcon::count","pcon::counter","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","pcon::dump","pcon::error","","","","","","","","","","","","","","","","","","","","pcon::serialize","","","","","","","","","","","","","","","","","","","","","","","","pcon::solid","","","","","","","","","","","","","","","pcon::spectrum","","","","","","","","","","","","","","","","","","","","","","",""],"d":["Command Line Interface declaration of project pcon","Run count command","Generic struct of counter and implementation for many type","Run dump command","Error struct of project pcon","Tools to serialize a Counter","Define Solid struct","Define Spectrum struct","Prompt COuNter, a short kmer counter.","SubCommand Count","Perform count of kmer","Output in csv mode","SubCommand Dump","Convert pcon native output in other format","Choose dump type","Output in bin mode","Output in solid mode","Enumeration of subcommand","Get abundance","Get abundance","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Returns the argument unchanged.","Returns the argument unchanged.","Returns the argument unchanged.","Returns the argument unchanged.","Returns the argument unchanged.","","","","","","","","","","","","","Get inputs","Get inputs","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","Get size of kmer","Get output","Get output","","Get quiet","Get record_buffer","SubCommand","Get timestamp granularity","","","","","","","","","","","","","","","","","","","","","","","","","","","Get verbosity level","Run count","A counter of kmer based on cocktail crate 2bit conversion, …","","","","","","Perform count on fasta input","Perform count on fasta input","Perform count on fasta input","Perform count on fasta input","Perform count on fasta input","","","","Returns the argument unchanged.","Create a new kmer by read a file","Create a new kmer by read a file","Create a new kmer by read a file","Create a new kmer by read a file","Create a new kmer by read a file","Get count of a kmer","Get count of a kmer","Get count of a kmer","Get count of a kmer","Get count of a kmer","Get count at on index","","Calls <code>U::from(self)</code>.","Get value of k","Create a new kmer Counter with kmer size equal to k","Create a new kmer Counter with kmer size equal to k","Create a new kmer Counter with kmer size equal to k","Create a new kmer Counter with kmer size equal to k","Create a new kmer Counter with kmer size equal to k","","Get raw data","Convert counter in serializer","","","","","Run dump","Error if we can’t convert a DumpTypeFromStr","Enum to manage error","Cost io error","Error in logging system configuration","Alias of result","Error durring loading count type not match","","","","","Returns the argument unchanged.","","","Calls <code>U::from(self)</code>.","","","","","","","Struct to serialize counter","","","Write kmer in csv format","Write kmer in csv format","Write kmer in csv format","Write kmer in csv format","Write kmer in csv format","Returns the argument unchanged.","Calls <code>U::from(self)</code>.","Create a new Serialize from a Counter","Write counter in pcon format","Write counter in pcon format","Write counter in pcon format","Write counter in pcon format","Write counter in pcon format","Convert counter in solid and write it","Convert counter in solid and write it","Convert counter in solid and write it","Convert counter in solid and write it","Convert counter in solid and write it","","","","A struct to store if a kmer is Solid or not. Only kmer …","","","Returns the argument unchanged.","Create a new Solid with count in <code>counter</code> only kmer upper …","Get the solidity status of <code>kmer</code>","Get the solidity status of a canonical <code>kmer</code>","Calls <code>U::from(self)</code>.","Get value of k","Create a new Solid for kmer size equal to <code>k</code>","Solidity status of <code>kmer</code> is set to <code>value</code>","Solidity status of a canonical<code>kmer</code> is set to <code>value</code>","","","","The first local minimum match with the intersection of …","Remove at least n percent of total kmer","Remove at most n percent of total kmer","More we remove kmer less we remove Erroneous kmer when …","A struct to represent kmer spectrum and usefull …","Based on Kmergenie we assume kmer spectrum is a mixture of …","","","","","","","Returns the argument unchanged.","Returns the argument unchanged.","Create a new Spectrum with count in <code>counter</code>","Found threshold matching with method","Calls <code>U::from(self)</code>.","Calls <code>U::from(self)</code>.","","","","","",""],"i":[0,0,0,0,0,0,0,0,0,0,11,5,0,11,0,5,5,0,1,3,8,1,3,8,1,3,11,11,8,11,5,1,3,8,11,5,1,3,5,5,5,8,8,5,8,11,5,1,3,8,11,5,1,3,8,11,1,3,8,11,1,3,8,1,3,11,3,1,8,11,5,1,3,1,1,3,5,8,1,8,8,5,5,8,11,5,1,3,8,11,5,1,3,8,11,5,1,3,8,11,1,3,8,11,1,3,5,8,0,0,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,0,39,0,39,39,0,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,0,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,38,0,45,45,45,45,45,45,45,45,45,45,45,45,45,45,46,46,46,46,0,0,47,46,47,46,46,46,47,46,47,47,47,46,47,46,47,46,47,46],"f":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,[1,2],[3,2],[4,4],[4,4],[4,4],[4,4],[4,4],[4,4],[4,4],[4,4],[[]],[[]],[[]],[[]],[[]],[[]],[[]],[[]],[[]],[[]],[5,5],[[]],[[5,5],6],[[],4],[[],4],[[5,5],7],[[8,9],10],[[11,9],10],[[5,9],10],[[1,9],10],[[3,9],10],[[]],[[]],[[]],[[]],[[]],[12,[[14,[8,13]]]],[12,[[14,[11,13]]]],[12,[[14,[1,13]]]],[12,[[14,[3,13]]]],[12,[[14,[8,13]]]],[12,[[14,[11,13]]]],[12,[[14,[1,13]]]],[12,[[14,[3,13]]]],[[],[[16,[15]]]],[[],[[16,[15]]]],[[],[[16,[15]]]],[17,7],[3,[[20,[[19,[18]]]]]],[1,[[20,[[19,[18]]]]]],[[]],[[]],[[]],[[]],[[]],[1,2],[1,21],[3,21],[[5,5],[[16,[6]]]],[8,7],[1,22],0,[8,23],[[]],[5,[[16,[24]]]],[[],14],[[],14],[[],14],[[],14],[[],14],[[],14],[[],14],[[],14],[[],14],[[],14],[[],25],[[],25],[[],25],[[],25],[[],25],[[8,12],[[14,[13]]]],[[11,12],[[14,[13]]]],[[1,12],[[14,[13]]]],[[3,12],[[14,[13]]]],[[8,12],[[14,[13]]]],[[11,12],[[14,[13]]]],[[1,12],[[14,[13]]]],[[3,12],[[14,[13]]]],[[]],[8,26],[1,20],0,[[]],[[]],[[[28,[27]]],[[28,[27]]]],[[]],[[[28,[29]],28],6],[[[28,[22]],[19,[18]],22]],[[[28,[30]],[19,[18]],22]],[[[28,[31]],[19,[18]],22]],[[[28,[32]],[19,[18]],22]],[[[28,[2]],[19,[18]],22]],[[],[[28,[33]]]],[[[28,[34]],28],7],[[[28,[35]],9],10],[[]],[[],[[20,[[28,[30]]]]]],[[],[[20,[[28,[22]]]]]],[[],[[20,[[28,[32]]]]]],[[],[[20,[[28,[2]]]]]],[[],[[20,[[28,[31]]]]]],[[[28,[32]],22],32],[[[28,[2]],22],2],[[[28,[22]],22],22],[[[28,[31]],22],31],[[[28,[30]],22],30],[[28,26]],[[[28,[36]]]],[[]],[28,2],[2,[[28,[32]]]],[2,[[28,[2]]]],[2,[[28,[22]]]],[2,[[28,[31]]]],[2,[[28,[30]]]],[[[28,[37]],28],[[16,[6]]]],[28],[28,38],[[]],[[],14],[[],14],[[],25],[3,20],0,0,0,0,0,0,[[]],[[]],[[39,9],10],[[39,9],10],[[]],[40,39],[41,39],[[]],[42],[39,[[16,[43]]]],[[],44],[[],14],[[],14],[[],25],0,[[]],[[]],[[[38,[32]],32],20],[[[38,[31]],31],20],[[[38,[30]],30],20],[[[38,[22]],22],20],[[[38,[2]],2],20],[[]],[[]],[28,38],[[[38,[30]]],20],[[[38,[22]]],20],[[[38,[31]]],20],[[[38,[32]]],20],[[[38,[2]]],20],[[[38,[22]],22],20],[[[38,[31]],31],20],[[[38,[32]],32],20],[[[38,[30]],30],20],[[[38,[2]],2],20],[[],14],[[],14],[[],25],0,[[]],[[]],[[]],[2,45],[[45,22],7],[[45,22],7],[[]],[45,2],[2,45],[[45,22,7]],[[45,22,7]],[[],14],[[],14],[[],25],0,0,0,0,0,0,[[]],[[]],[[]],[[]],[[46,46],7],[[46,9],10],[[]],[[]],[[],47],[[47,46,48],[[16,[2]]]],[[]],[[]],[[],14],[[],14],[[],14],[[],14],[[],25],[[],25]],"p":[[3,"Count"],[15,"u8"],[3,"Dump"],[3,"Command"],[4,"DumpType"],[4,"Ordering"],[15,"bool"],[3,"Command"],[3,"Formatter"],[6,"Result"],[4,"SubCommand"],[3,"ArgMatches"],[6,"Error"],[4,"Result"],[3,"Id"],[4,"Option"],[15,"str"],[8,"BufRead"],[3,"Box"],[6,"Result"],[3,"Vec"],[15,"u64"],[4,"Timestamp"],[3,"PossibleValue"],[3,"TypeId"],[15,"usize"],[8,"Clone"],[3,"Counter"],[8,"Ord"],[15,"u32"],[15,"u128"],[15,"u16"],[8,"Default"],[8,"PartialEq"],[8,"Debug"],[8,"Hash"],[8,"PartialOrd"],[3,"Serialize"],[4,"Error"],[3,"SetLoggerError"],[3,"Error"],[3,"Demand"],[8,"Error"],[3,"String"],[3,"Solid"],[4,"ThresholdMethod"],[3,"Spectrum"],[15,"f64"]]}\
}');
if (typeof window !== 'undefined' && window.initSearch) {window.initSearch(searchIndex)};
if (typeof exports !== 'undefined') {exports.searchIndex = searchIndex};