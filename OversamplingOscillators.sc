SawOS : UGen {
	*ar { |freq=440, phase=0, oversample=1, mul = 1, add = 0|
		^this.multiNew('audio', freq, phase, oversample).madd(mul, add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

VarSawOS : UGen {
	*ar { |freq=440, iphase=0, width = 0.5, oversample=1, mul = 1, add = 0|
		^this.multiNew('audio', freq, iphase, width, oversample).madd(mul, add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

SawPn : UGen {
	*ar { |freq=440, phase=0, mul=1, add=0|
		^this.multiNew('audio', freq, phase).madd(mul,add);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}