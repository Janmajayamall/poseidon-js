const assert = require("assert");
const { getCurveFromName, utils } = require("ffjavascript");
const fs = require("fs");

class Poseidon {
	constructor(specPath) {
		this.specPath = specPath;

		this.spec = undefined;
		this.F = undefined;
		this.curve = undefined;

		this.state = [];
		this.absorbing = [];
	}

	async parseSpec() {
		if (this.F == undefined) throw new Error("Curve not loaded!");
		const F = this.F;

		this.spec = await new Promise((resolve, reject) => {
			fs.readFile(this.specPath, "utf-8", (error, data) => {
				if (error) {
					reject(`spec parsing failed with ${error}`);
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
						roundConstants.push(F.fromRprLE(roundConstantsRaw[i]));
					}
					constants.start.push(roundConstants);
				}
				constants.end = [];
				for (let j = 0; j < raw.constants.end.length; j++) {
					let roundConstantsRaw = raw.constants.end[j];
					let roundConstants = [];
					for (let i = 0; i < roundConstantsRaw.length; i++) {
						roundConstants.push(F.fromRprLE(roundConstantsRaw[i]));
					}
					constants.end.push(roundConstants);
				}
				constants.partial = [];
				for (let j = 0; j < raw.constants.partial.length; j++) {
					constants.partial.push(
						F.fromRprLE(raw.constants.partial[j])
					);
				}

				let mdsMatrices = {};
				mdsMatrices.mds = [];
				for (let j = 0; j < raw.mdsMatrices.mds.length; j++) {
					let rowRaw = raw.mdsMatrices.mds[j];
					let row = [];
					for (let i = 0; i < rowRaw.length; i++) {
						row.push(F.fromRprLE(rowRaw[i]));
					}
					mdsMatrices.mds.push(row);
				}
				mdsMatrices.preSparseMds = [];
				for (let j = 0; j < raw.mdsMatrices.preSparseMds.length; j++) {
					let rowRaw = raw.mdsMatrices.preSparseMds[j];
					let row = [];
					for (let i = 0; i < rowRaw.length; i++) {
						row.push(F.fromRprLE(rowRaw[i]));
					}
					mdsMatrices.preSparseMds.push(row);
				}
				mdsMatrices.sparseMatrices = [];
				for (
					let j = 0;
					j < raw.mdsMatrices.sparseMatrices.length;
					j++
				) {
					let sparseMatrixRaw = raw.mdsMatrices.sparseMatrices[j];
					let sparseMatrix = {};
					sparseMatrix.row = [];
					for (let i = 0; i < sparseMatrixRaw.row.length; i++) {
						sparseMatrix.row.push(
							F.fromRprLE(sparseMatrixRaw.row[i])
						);
					}
					sparseMatrix.colHat = [];
					for (let i = 0; i < sparseMatrixRaw.colHat.length; i++) {
						sparseMatrix.colHat.push(
							F.fromRprLE(sparseMatrixRaw.colHat[i])
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

	async loadCurve() {
		this.curve = await getCurveFromName("bn128", true);
		this.F = this.curve.Fr;
	}

	loadState() {
		if (this.F == undefined) throw new Error("Curve not loaded");
		if (this.spec == undefined) throw new Error("Spec not loaded");

		const F = this.F;

		let state = [];
		for (let i = 0; i < this.spec.t; i++) {
			state.push(F.e(0));
		}
		// 2 ** 64
		state[0] = F.e(18446744073709551616);

		this.state = state;
	}

	permute() {
		const spec = this.spec;
		const F = this.F;
		let state = this.state;

		// helpers
		const sBoxFull = (state) => {
			return state.map((val) => {
				return F.mul(val, F.square(F.square(val, val)));
			});
		};
		const matrixMulVector = (matrix, vector) => {
			let res = [];
			for (let i = 0; i < matrix.length; i++) {
				res.push(
					vector.reduce((acc, v, j) => {
						return F.add(acc, F.mul(v, matrix[i][j]));
					}, F.e(0))
				);
			}
			return res;
		};
		const sBoxPartial = (state) => {
			state[0] = F.mul(state[0], F.square(F.square(state[0], state[0])));
			return state;
		};
		const addRoundConstants = (state, constants) => {
			return state.map((val, i) => {
				return F.add(val, constants[i]);
			});
		};

		// T should always equal state length
		if (state.length != spec.t) throw new Error("T != state.length");

		let rF = spec.rF / 2;

		// First half full rounds
		state = addRoundConstants(state, spec.constants.start[0]);
		for (let i = 1; i < rF; i++) {
			state = sBoxFull(state);
			state = addRoundConstants(state, spec.constants.start[i]);
			state = matrixMulVector(spec.mdsMatrices.mds, state);
		}
		state = sBoxFull(state);
		state = addRoundConstants(
			state,
			spec.constants.start[spec.constants.start.length - 1]
		);
		state = matrixMulVector(spec.mdsMatrices.preSparseMds, state);

		// Partial rounds
		for (let i = 0; i < spec.constants.partial.length; i++) {
			state = sBoxPartial(state);
			state[0] = F.add(state[0], spec.constants.partial[i]);

			// apply sparse matrix mds to state
			let sparseMatrix = spec.mdsMatrices.sparseMatrices[i];
			let tmp = state[0];
			state[0] = sparseMatrix.row.reduce((acc, v, index) => {
				return F.add(acc, F.mul(v, state[index]));
			}, F.e(0));
			for (let j = 1; j < state.length; j++) {
				state[j] = F.add(
					F.mul(sparseMatrix.colHat[j - 1], tmp),
					state[j]
				);
			}
		}

		// Second half full rounds
		for (let i = 0; i < spec.constants.end.length; i++) {
			state = sBoxFull(state);
			state = addRoundConstants(state, spec.constants.end[i]);
			state = matrixMulVector(spec.mdsMatrices.mds, state);
		}
		state = sBoxFull(state);
		state = matrixMulVector(spec.mdsMatrices.mds, state);

		this.state = state;
	}

	updateFromRprLE(inputs) {
		const F = this.F;

		if (!Array.isArray(inputs)) throw new Error("Input must be an array!");

		inputs = inputs.map((input) => {
			assert(
				(Array.isArray(input) && input.length == F.n8) ||
					input.byteLength == F.n8
			);

			return F.fromRprLE(input);
		});

		this.update(inputs);
	}

	update(inputs) {
		const F = this.F;

		assert(
			inputs.reduce((acc, v, i) => {
				return v.byteLength == F.n8 && acc;
			}),
			true
		);

		// add inputs absorption line
		let absorbing = this.absorbing.concat(inputs);

		let rate = this.spec.rate;
		let pInputs = [];
		let k = 0;
		for (let i = 0; i < absorbing.length; i++) {
			pInputs.push(absorbing[i]);

			if ((i != 0 || rate == 1) && (i + 1) % rate == 0) {
				// add pInputs to state
				assert(this.state.length - 1 == pInputs.length);
				this.state = this.state.map((v, i) => {
					if (i != 0) {
						return F.add(v, pInputs[i - 1]);
					}
					return v;
				});

				// permute
				this.permute();

				pInputs = [];
				k = i;
			}
		}

		this.absorbing = absorbing.slice(k + 1);
		assert(this.absorbing.length < rate);
	}

	squeeze() {
		const F = this.F;

		let rate = this.spec.rate;
		let absorbing = this.absorbing;
		assert(absorbing.length < rate);

		absorbing.push(F.e(1));

		for (let i = 0; i < absorbing.length; i++) {
			this.state[i + 1] = F.add(this.state[i + 1], absorbing[i]);
		}

		this.permute();

		this.absorbing = [];

		return this.state[1];
	}

	printFormatted(value, text) {
		text = text || "";
		const F = this.F;
		// let tmp = new Uint8Array(F.n8);
		// F.toRprLE(tmp, 0, value);
		console.log(F.toString(value, 16), text);
	}

	printState(state) {
		state.forEach((s) => {
			this.printFormatted(s);
		});
	}
}

async function main() {
	let poseidon = new Poseidon("./tmp/spec_values.json");
	await poseidon.loadCurve();
	await poseidon.parseSpec();
	poseidon.loadState();

	let input = [
		[
			31, 74, 114, 51, 235, 141, 133, 7, 188, 101, 67, 146, 252, 94, 5,
			219, 173, 155, 238, 219, 228, 155, 159, 93, 44, 131, 193, 94, 141,
			243, 137, 37,
		],
		[
			31, 74, 114, 51, 235, 141, 133, 7, 188, 101, 67, 146, 252, 94, 5,
			219, 173, 155, 238, 219, 228, 155, 159, 93, 44, 131, 193, 94, 141,
			243, 137, 37,
		],
	];

	poseidon.updateFromRprLE(input);
	let hash = poseidon.squeeze();
	poseidon.printFormatted(hash, " final hash");
}

main()
	.then(() => {})
	.catch((e) => {
		console.log(e);
	});
