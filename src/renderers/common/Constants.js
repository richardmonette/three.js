export const AttributeType = {
	VERTEX: 1,
	INDEX: 2,
	STORAGE: 3,
	INDIRECT: 4
};

// size of a chunk in bytes (STD140 layout)

export const GPU_CHUNK_BYTES = 16;

// max gap in elements between two buffer update ranges that are still uploaded with a single call

export const UPDATE_RANGE_MERGE_GAP = 64;
