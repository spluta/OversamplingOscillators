SawOS : UGen {
	*ar { |freq=440, phase=0, oversample=1, mul = 1, add = 0|
		if(freq.rate!='audio'){freq = K2A.ar(freq)};
		if(phase.rate!='audio'){phase = K2A.ar(phase)};
		^this.multiNew('audio', freq, phase, oversample).madd(mul, add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

SinOscOS : UGen {
	*ar { |freq=440, phase=0, oversample=1, mul = 1, add = 0|
		if(freq.rate!='audio'){freq = K2A.ar(freq)};
		if(phase.rate!='audio'){phase = K2A.ar(phase)};
		^this.multiNew('audio', freq, phase, oversample).madd(mul, add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

PMOscOS : UGen {//CarFreq, ModFreq, PMMul, PMModPhase, OverSample
	*ar { |carFreq=440, modFreq=220, pmMul=1, pmModPhase=0, oversample=1, mul = 1, add = 0|
		if(carFreq.rate!='audio'){carFreq = K2A.ar(carFreq)};
		if(modFreq.rate!='audio'){modFreq = K2A.ar(modFreq)};
		if(pmMul.rate!='audio'){pmMul = K2A.ar(pmMul)};
		if(pmModPhase.rate!='audio'){pmModPhase = K2A.ar(pmModPhase)};
		^this.multiNew('audio', carFreq, modFreq, pmMul, pmModPhase, oversample).madd(mul, add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

FM7OS : MultiOutUGen  {
	*numOperators { ^6 }
	*ar { | ctlMatrix, modMatrix, oversample=4, mul=1, add=0 |
		ctlMatrix = (ctlMatrix ?? { this.controlMatrix })
		.collect{|item|
			if(item[0].rate!='audio'){item[0] = K2A.ar(item[0])};
			if(item[2].rate!='audio'){item[2] = K2A.ar(item[2])};
			item
		}.flatten(1);
		modMatrix = (modMatrix ?? { this.modMatrix }).flatten(1)
		.collect{|item|
			if(item.rate!='audio'){item = K2A.ar(item)};
			item
		};
		^this.multiNewList(
			['audio']
			++ ctlMatrix
			++ modMatrix
			++ oversample
		).madd(mul,add)
	}

	init { | ... args |
		inputs = args;
		^this.initOutputs(this.class.numOperators, rate)
	}
	checkInputs {
		^this.checkValidInputs;
	}
}

FM7aOS : MultiOutUGen  {
	*ar { | ctlMatrix, modMatrix, synthTypes, oversample=4, mul=1, add=0 |
		ctlMatrix = (ctlMatrix ?? { this.controlMatrix })
		.collect{|item|
			if(item.rate!='audio'){item = K2A.ar(item)};
			item
		};
		synthTypes.do{|type,i| modMatrix[i][i] = type};
		modMatrix.postln;
		modMatrix = modMatrix.flatten(1)
		.collect{|item|
			if(item.rate!='audio'){item = K2A.ar(item)};
			item
		};
		^this.multiNewList(
			['audio']
			++ ctlMatrix
			++ modMatrix
			++ oversample
		).madd(mul,add)
	}

	init { | ... args |
		inputs = args;
		^this.initOutputs(4, rate)
	}
	checkInputs {
		^this.checkValidInputs;
	}
}

FM7bOS : MultiOutUGen  {
*ar { | ctlMatrix, modMatrix, synthTypes, oversample=4, mul=1, add=0 |
		ctlMatrix = (ctlMatrix ?? { this.controlMatrix })
		.collect{|item|
			if(item.rate!='audio'){item = K2A.ar(item)};
			item
		};
		synthTypes.do{|type,i| modMatrix[i][i] = type};
		modMatrix.postln;
		modMatrix = modMatrix.flatten(1)
		.collect{|item|
			if(item.rate!='audio'){item = K2A.ar(item)};
			item
		};
		^this.multiNewList(
			['audio']
			++ ctlMatrix
			++ modMatrix
			++ oversample
		).madd(mul,add)
	}

	init { | ... args |
		inputs = args;
		^this.initOutputs(4, rate)
	}
	checkInputs {
		^this.checkValidInputs;
	}
}

PM7OS : MultiOutUGen {
	*numOperators { ^6 }
	*numControls { ^3 }
	*controlArraySize { ^this.numControls * this.numOperators }
	*modArraySize { ^this.numOperators.squared }

	*controlMatrix { | ... args |
		var matrix;
		matrix = Array.fill2D(this.numOperators, this.numControls, 0);
		args.do { | x |
			matrix[x[0]] = x.copyToEnd(1);
		};
		^matrix
	}
	*modMatrix { | ... args |
		var matrix;
		matrix = Array.fill2D(this.numOperators, this.numOperators, 0);
		args.do { | x |
			matrix[x[0]][x[1]] = x[2];
		};
		^matrix
	}

	*algoSpec { | algo, feedback=0.0 |
		var matrix;
		matrix = this.modMatrix;
		^[
			{[
				this.modMatrix(
					[0, 1, 1],
					[2, 3, 1],
					[3, 4, 1],
					[4, 5, 1],
					[5, 5, feedback]
				),  [0, 2]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[1, 1, feedback],
					[2, 3, 1],
					[3, 4, 1],
					[4, 5, 1]
				),  [0, 2]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[1, 2, 1],
					[3, 4, 1],
					[4, 5, 1],
					[5, 5, feedback]
				),  [0, 3]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[1, 2, 1],
					[3, 4, 1],
					[4, 5, 1],
					[5, 3, feedback]
				),  [0, 3]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[2, 3, 1],
					[5, 5, feedback],
				),  [0, 2, 4]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[2, 3, 1],
					[5, 4, feedback],
				),  [0, 2, 4]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[2, 3, 1],
					[2, 4, 1],
					[4, 5, 1],
					[5, 5, feedback],
				),  [0, 2]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[2, 3, 1],
					[2, 4, 1],
					[4, 5, 1],
					[3, 3, feedback],
				),  [0, 2]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[2, 3, 1],
					[2, 4, 1],
					[4, 5, 1],
					[1, 1, feedback],
				),  [0, 2]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[1, 2, 1],
					[3, 4, 1],
					[3, 5, 1],
					[2, 2, feedback],
				),  [0, 3]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[1, 2, 1],
					[3, 4, 1],
					[3, 5, 1],
					[5, 5, feedback],
				),  [0, 3]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[1, 2, 1],
					[2, 3, 1],
					[2, 4, 1],
					[2, 5, 1],
					[2, 2, feedback],
				),  [0, 2]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[1, 2, 1],
					[2, 3, 1],
					[2, 4, 1],
					[2, 5, 1],
					[5, 5, feedback],
				),  [0, 2]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[2, 3, 1],
					[2, 3, 1],
					[3, 4, 1],
					[3, 5, 1],
					[5, 5, feedback],
				),  [0, 2]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[2, 3, 1],
					[2, 3, 1],
					[3, 4, 1],
					[3, 5, 1],
					[1, 1, feedback],
				),  [0, 2]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[0, 2, 1],
					[0, 4, 1],
					[2, 3, 1],
					[4, 5, 1],
					[5, 5, feedback],
				),  [0]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[0, 2, 1],
					[0, 4, 1],
					[2, 3, 1],
					[4, 5, 1],
					[1, 1, feedback],
				),  [0]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[0, 2, 1],
					[0, 3, 1],
					[3, 4, 1],
					[4, 5, 1],
					[2, 2, feedback],
				),  [0]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[1, 2, 1],
					[3, 5, 1],
					[4, 5, 1],
					[5, 5, feedback],
				),  [0, 3, 4]
			]},
			{[
				this.modMatrix(
					[0, 2, 1],
					[1, 2, 1],
					[3, 4, 1],
					[3, 5, 1],
					[2, 2, feedback],
				),  [0, 1, 3]
			]},
			{[
				this.modMatrix(
					[0, 2, 1],
					[1, 2, 1],
					[3, 5, 1],
					[4, 5, 1],
					[2, 2, feedback],
				),  [0, 1, 3, 4]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[2, 5, 1],
					[3, 5, 1],
					[4, 5, 1],
					[5, 5, feedback],
				),  [0, 2, 3, 4]
			]},
			{[
				this.modMatrix(
					[1, 2, 1],
					[3, 5, 1],
					[4, 5, 1],
					[5, 5, feedback],
				),  [0, 1, 3, 4]
			]},
			{[
				this.modMatrix(
					[2, 5, 1],
					[3, 5, 1],
					[4, 5, 1],
					[5, 5, feedback],
				),  [0, 1, 2, 3, 4]
			]},
			{[
				this.modMatrix(
					[3, 5, 1],
					[4, 5, 1],
					[5, 5, feedback],
				),  [0, 1, 2, 3, 4]
			]},
			{[
				this.modMatrix(
					[1, 2, 1],
					[3, 4, 1],
					[3, 5, 1],
					[5, 5, feedback],
				),  [0, 1, 3]
			]},
			{[
				this.modMatrix(
					[1, 2, 1],
					[3, 4, 1],
					[3, 5, 1],
					[2, 2, feedback],
				),  [0, 1, 3]
			]},
			{[
				this.modMatrix(
					[0, 1, 1],
					[2, 3, 1],
					[3, 4, 1],
					[4, 4, feedback],
				),  [0, 2, 5]
			]},
			{[
				this.modMatrix(
					[2, 3, 1],
					[4, 5, 1],
					[5, 5, feedback],
				),  [0, 1, 2, 4]
			]},
			{[
				this.modMatrix(
					[2, 3, 1],
					[3, 4, 1],
					[4, 4, feedback],
				),  [0, 1, 2, 5]
			]},
			{[
				this.modMatrix(
					[4, 5, 1],
					[5, 5, feedback],
				),  [0, 1, 2, 3, 4]
			]},
			{[
				this.modMatrix(
					[5, 5, feedback],
				),  [0, 1, 2, 3, 4, 5]
			]},
		].at(algo).value*[1/pi, 1]  //key is that it is divided by pi
	}

	*ar { | ctlMatrix, modMatrix, oversample=4, mul=1, add=0 |
		ctlMatrix = (ctlMatrix ?? { this.controlMatrix })
		.collect{|item|
			if(item[0].rate!='audio'){item[0] = K2A.ar(item[0])};
			if(item[2].rate!='audio'){item[2] = K2A.ar(item[2])};
			item
		}.flatten(1);
		modMatrix = (modMatrix ?? { this.modMatrix }).flatten(1)
		.collect{|item|
			if(item.rate!='audio'){item = K2A.ar(item)};
			item
		};
		^this.multiNewList(
			['audio']
			++ ctlMatrix
			++ modMatrix
			++ oversample
		).madd(mul,add)
	}

	*arAlgo { | algo=0, ctlMatrix, feedback=0.0, oversample=4 |
		var modMatrix, channels;
		#modMatrix, channels = this.algoSpec(algo, feedback);
		^this.ar(ctlMatrix, modMatrix, oversample).slice(channels)
	}

	init { | ... args |
		inputs = args;
		^this.initOutputs(this.class.numOperators, rate)
	}
	checkInputs {
		^this.checkValidInputs;
	}
}


SquareOS : UGen {
	*ar { |freq=440, phase=0, width = 0.5, oversample=1, mul = 1, add = 0|
		if(freq.rate!='audio'){freq = K2A.ar(freq)};
		if(phase.rate!='audio'){phase = K2A.ar(phase)};
		if(width.rate!='audio'){width = K2A.ar(width)};
		^this.multiNew('audio', freq, phase, width, oversample).madd(mul, add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

TriOS : UGen {
	*ar { |freq=440, phase=0, oversample=1, mul = 1, add = 0|
		if(freq.rate!='audio'){freq = K2A.ar(freq)};
		if(phase.rate!='audio'){phase = K2A.ar(phase)};
		^this.multiNew('audio', freq, phase, oversample).madd(mul, add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

VarSawOS : UGen {
	*ar { |freq=440, phase=0, width = 0.5, oversample=1, mul = 1, add = 0|
		if(freq.rate!='audio'){freq = K2A.ar(freq)};
		if(width.rate!='audio'){width = K2A.ar(width)};
		if(phase.rate!='audio'){phase = K2A.ar(phase)};
		^this.multiNew('audio', freq, phase, width, oversample).madd(mul, add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

BuchlaFoldOS : UGen {
	*ar { |sig, amp = 1, oversample=1, mul = 1, add = 0|
		if(sig.rate!='audio'){sig = K2A.ar(sig)};
		if(amp.rate!='audio'){amp = K2A.ar(amp)};
		^this.multiNew('audio', sig, amp, oversample).madd(mul, add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

SergeFoldOS : UGen {
	*ar { |sig, amp = 1, oversample=1, mul = 1, add = 0|
		if(sig.rate!='audio'){sig = K2A.ar(sig)};
		if(amp.rate!='audio'){amp = K2A.ar(amp)};
		^this.multiNew('audio', sig, amp, oversample).madd(mul, add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

SawBL : UGen {
	*ar { |freq=440, iphase=0, mul=1, add=0|
		if(freq.rate!='audio'){freq = K2A.ar(freq)};
		^this.multiNew('audio', freq, iphase).madd(mul,add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

SquareBL : UGen {
	*ar { |freq=440, phase=0, width=0.5, mul=1, add=0|
		if(freq.rate!='audio'){freq = K2A.ar(freq)};
		if(width.rate!='audio'){width = K2A.ar(width)};
		^this.multiNew('audio', freq, phase, width).madd(mul,add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

TriBL : UGen {
	*ar { |freq=440, phase=0, mul=1, add=0|
		var width = 0.5;
		if(freq.rate!='audio'){freq = K2A.ar(freq)};
		if(width.rate!='audio'){width = K2A.ar(width)};
		^this.multiNew('audio', freq, phase, width).madd(mul,add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

ImpulseBL : UGen {
	*ar { |freq=440, mul=1, add=0|
		if(freq.rate!='audio'){freq = K2A.ar(freq)};
		^this.multiNew('audio', freq).madd(mul,add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}