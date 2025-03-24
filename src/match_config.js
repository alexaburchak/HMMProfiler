import * as node_fs from "node:fs";
import yargs from "yargs/yargs";

/**
 * Represents a pair of input files for the pipeline
 * @typedef {object} input_list
 * @property {string} query_path - Path to the query sequences (can be a FASTA file or a single sequence string)
 * @property {string} model_path - Path to hmm for trimming query sequences
 * @property {string} csv_path - Path to counts csv file to be searched
 * @property {string} output_path - Path to output csv of identified matches
 */

/**
 * @typedef {object} PipelineConfig
 * @property {number} max_LD - Maximum levenshtein distance
 * @property {input_list[]} input_list - List of query sequences to process and count csv files for searching
 */

/**
 * Validates whether the given input is an input_list object
 * @param {any} maybe_input_list
 * @returns {boolean}
 */
function is_input(maybe_input_list) {
	if (typeof maybe_input_list !== "object") {
		return false;
	}
	if (maybe_input_list === null) {
		return false;
	}
	const query_path = maybe_input_list.query_path;
	if (typeof query_path !== "string") {
		return false;
	}
	const model_path = maybe_input_list.model_path;
	if (typeof model_path !== "string") {
		return false;
	}
	const csv_path = maybe_input_list.csv_path;
	if (typeof csv_path !== "string") {
		return false;
	}
	const output_path = maybe_input_list.output_path;
	if (typeof output_path !== "string") {
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
	const max_LD = maybe_config.max_LD;
	if (typeof max_LD !== "number") {
		return false;
	}
	const input_list = maybe_config.input_list;
	if (!Array.isArray(input_list)) {
		return false;
	}
	for (const input of input_list) {
		if (!is_input(input)) {
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
