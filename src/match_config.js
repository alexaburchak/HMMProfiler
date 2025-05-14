import * as node_fs from "node:fs";
import yargs from "yargs/yargs";

/**
 * @typedef {object} MatchesConfig
 * @property {{name: string, sequences: string[]}[]} queryEntries - Array of query objects with names (string) and sequences (array of strings).
 * @property {number} max_LD - Maximum Levenshtein distance for matching sequences.
 * @property {number} hmm_coverage - Minimum required model coverage for hmmsearch hits.
 * @property {{name: string, model_paths: string[], counts_path: string}[]} libraries -  Array of objects with library name (string), model_paths (array of strings), and counts_path (string).
 * @property {string} output_path - Path to write CSV output of all detected matches.
 */

/**
 * Validates whether the given input is a PipelineConfig object
 * @param {any} maybe_config
 * @returns {boolean}
 */
function is_config(maybe_config) {
	if (typeof maybe_config !== "object") {
		return false;
	}
	if (maybe_config === null) {
		return false;
	}

	// Validate queryEntries
	if (
		!Array.isArray(maybe_config.queryEntries) ||
		maybe_config.queryEntries.length === 0
	) {
		return false;
	}

	for (const entry of maybe_config.queryEntries) {
		if (typeof entry !== "object" || entry === null) {
			return false;
		}
		const { name, sequences } = entry;

		if (typeof name !== "string") {
			return false;
		}
		if (
			!Array.isArray(sequences) ||
			sequences.length === 0 ||
			!sequences.every((seq) => typeof seq === "string")
		) {
			return false;
		}
	}

	// Validate libraries
	if (
		!Array.isArray(maybe_config.libraries) ||
		maybe_config.libraries.length === 0
	) {
		return false;
	}

	for (const library of maybe_config.libraries) {
		if (typeof library !== "object" || library === null) {
			return false;
		}
		const { name, model_paths, counts_path } = library;

		if (typeof name !== "string") {
			return false;
		}
		if (
			!Array.isArray(model_paths) ||
			model_paths.length === 0 ||
			!model_paths.every((path) => typeof path === "string")
		) {
			return false;
		}
		if (typeof counts_path !== "string") {
			return false;
		}
	}

	// Validate output path
	const output_path = maybe_config.output_path;
	if (typeof output_path !== "string") {
		return false;
	}

	// Validate levenshtein distance filter
	const max_LD = maybe_config.max_LD;
	if (typeof max_LD !== "number") {
		return false;
	}

	// Validate hmm coverage filter 
	const hmm_coverage = maybe_config.hmm_coverage;
	if (typeof hmm_coverage !== "number") {
		return false;
	}
	return true;
}

/**
 * Reads + parses a JSON configuration file from given path
 * @param {string} config_path
 * @returns {MatchesConfig | null}
 */
function get_config_by_path(config_path) {
	const config_buffer = node_fs.readFileSync(config_path);
	const config_string = config_buffer.toString();
	try {
		const config = JSON.parse(config_string);
		if (!is_config(config)) {
			return null;
		}
		return config;
	} catch (error) {
		return null;
	}
}

/**
 * Parses command-line arguments to extract config file path
 * @returns {Promise<string>}
 */
async function get_config_path_by_args() {
	const argv = await yargs(process.argv.slice(2))
		.usage("Usage: $0 [options]")
		.option("c", {
			describe: "Config path",
			type: "string",
			demandOption: true,
		})
		.help("h")
		.alias("h", "help").argv;

	const config_path = argv.c;
	return config_path;
}

export { get_config_path_by_args, get_config_by_path };
