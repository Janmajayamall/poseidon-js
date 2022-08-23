const { getCurveFromName } = require("ffjavascript");
const fs = require("fs");

async function parseSpec(F) {
	return new Promise((resolve, reject) => {
		fs.readFile("./spec_values.json", "utf-8", (error, data) => {
			if (error) {
				reject(error);
			}

			let raw = JSON.parse(data);

			// TODO: we need to check that length of `start` and `end` is correct.
			// That is to say that they sum upto `rF` and the division of full rounds
			// between them is correct. This is important since the implementation
			// relies on these arrays for calculating full rounds (for ex, when performing
			// second half of full rounds we simply iterate over `end` array of round constants).
			// Also make sure that length of `partial` matches expected number of partial
			// rounds.
			let constants = {};
			constants.start = [];
			for (let j = 0; j < raw.constants.start.length; j++) {
				let roundConstantsRaw = raw.constants.start[j];
				let roundConstants = [];
				for (let i = 0; i < roundConstantsRaw.length; i++) {
					roundConstants.push(
						F.e(new Uint8Array(roundConstantsRaw[i]))
					);
				}
				constants.start.push(roundConstants);
			}
			constants.end = [];
			for (let j = 0; j < raw.constants.end.length; j++) {
				let roundConstantsRaw = raw.constants.end[j];
				let roundConstants = [];
				for (let i = 0; i < roundConstantsRaw.length; i++) {
					roundConstants.push(
						F.e(new Uint8Array(roundConstantsRaw[i]))
					);
				}
				constants.end.push(roundConstants);
			}
			constants.partial = [];
			for (let j = 0; j < raw.constants.partial.length; j++) {
				constants.partial.push(
					F.e(new Uint8Array(raw.constants.partial[j]))
				);
			}

			let mdsMatrices = {};
			mdsMatrices.mds = [];
			for (let j = 0; j < raw.mdsMatrices.mds.length; j++) {
				let rowRaw = raw.mdsMatrices.mds[j];
				let row = [];
				for (let i = 0; i < rowRaw.length; i++) {
					row.push(F.e(new Uint8Array(rowRaw[i])));
				}
				mdsMatrices.mds.push(row);
			}
			mdsMatrices.preSparseMds = [];
			for (let j = 0; j < raw.mdsMatrices.preSparseMds.length; j++) {
				let rowRaw = raw.mdsMatrices.preSparseMds[j];
				let row = [];
				for (let i = 0; i < rowRaw.length; i++) {
					row.push(F.e(new Uint8Array(rowRaw[i])));
				}
				mdsMatrices.preSparseMds.push(row);
			}
			mdsMatrices.sparseMatrices = [];
			for (let j = 0; j < raw.mdsMatrices.sparseMatrices.length; j++) {
				let sparseMatrixRaw = raw.mdsMatrices.sparseMatrices[j];
				let sparseMatrix = {};
				sparseMatrix.row = [];
				for (let i = 0; i < sparseMatrixRaw.row.length; i++) {
					sparseMatrix.row.push(
						F.e(new Uint8Array(sparseMatrixRaw.row[i]))
					);
				}
				sparseMatrix.colHat = [];
				for (let i = 0; i < sparseMatrixRaw.colHat.length; i++) {
					sparseMatrix.colHat.push(
						F.e(new Uint8Array(sparseMatrixRaw.colHat[i]))
					);
				}
				mdsMatrices.sparseMatrices.push(sparseMatrix);
			}

			resolve({
				constants,
				mdsMatrices,
				rF: raw.rF,
				rP: raw.rP,
				t: raw.t,
				rate: raw.rate,
			});
		});
	});
}

function permute(state, spec) {
	const sBoxFull = (state) => {
		return state.map((val) => {
			return F.mul(val, F.square(F.square(val, val)));
		});
	};
	const addRoundConstants = (state, constants) => {
		return state.map((val, i) => {
			F.add(val, constants[i]);
		});
	};

	let rF = spec.rF / 2;

	// First half full rounds
	state = addRoundConstants(state, spec.constants.start[0]);
	for (let i = 1; i < rF - 1; i++) {
		state = sBoxFull(state);
		state = addRoundConstants(state, spec.constants.start[i]);
		// apply MDS
	}
	state = sBoxFull(state);
	state = addRoundConstants(
		state,
		spec.constants.start[spec.constants.start.length - 1]
	);
	// apple MDS

	// Partial rounds

	// Second half full rounds
	for (let i = 0; i < spec.constants.end.length; i++) {
		state = sBoxFull(state);
		state = addRoundConstants(state, spec.constants.end[i]);
		// apply MDS
	}
	state = sBoxFull(state);
	// apply MDS
}

async function poseidon() {
	const bn128 = await getCurveFromName("bn128", true);
	const F = bn128.Fr;
	let spec = await parseSpec(F);
}

poseidon()
	.then(() => {})
	.catch((e) => {
		console.log(e);
	});
