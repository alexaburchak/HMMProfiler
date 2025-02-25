import * as node_fs from "node:fs";
import yargs from "yargs/yargs";

/**
 * Represents a pair of input files for the pipeline
 * @typedef {object} InputPair
 * @property {string} fastq_path - Path to the FASTQ file
 * @property {string} model_path - Path to the corresponding HMM model
 */

/**
 * @typedef {object} PipelineConfig
 * @property {string} matches_outpath - Path for hmmer matches output file
 * @property {string} counts_outpath - Path for VH/VL pair counts output file 
 * @property {number} min_quality - Minimum quality threshold for filtering reads
 * @property {InputPair[]} input_pairs - List of FASTQ and model file pairs to process
 */

/**
 * Validates whether the given input is an InputPair object
 * @param {any} maybe_input_pair
 * @returns {boolean}
 */
function is_input_pair(maybe_input_pair) {
	if (typeof maybe_input_pair !== "object") {
		return false;
	}
	if (maybe_input_pair === null) {
		return false;
	}
	const fastq_path = maybe_input_pair.fastq_path;
	if (typeof fastq_path !== "string") {
		return false;
	}
	const model_path = maybe_input_pair.model_path;
	if (typeof model_path !== "string") {
		return false;
	}

	return true;
}

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
	const min_quality = maybe_config.min_quality;
	if (typeof min_quality !== "number") {
		return false;
	}
	const matches_outpath = maybe_config.matches_outpath;
	if (typeof matches_outpath !== "string") {
		return false;
	}
	const counts_outpath = maybe_config.counts_outpath;
	if (typeof counts_outpath !== "string") {
		return false;
	}
	const input_pairs = maybe_config.input_pairs;
	if (!Array.isArray(input_pairs)) {
		return false;
	}
	for (const input_pair of input_pairs) {
		if (!is_input_pair(input_pair)) {
			return false;
		}
	}

	return true;
}

/**
 * Reads + parses a JSON configuration file from given path
 * @param {string} config_path
 * @returns {PipelineConfig | null}
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
