SawOS : UGen {
	*ar { |freq=440, phase=0, oversample=1|
		^this.multiNew('audio', freq, phase, oversample);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

VarSawOS : UGen {
	*ar { |freq=440, phase=0, width = 0.5, oversample=1|
		^this.multiNew('audio', freq, phase, width, oversample);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}

SawPn : UGen {
	*ar { |freq=440, phase=0, oversample=1|
		^this.multiNew('audio', freq, phase, oversample);
	}

	checkInputs {
		/* TODO */
		^this.checkValidInputs;
	}
}