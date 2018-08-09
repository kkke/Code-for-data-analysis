/*
 * Macro template to process multiple images in a folder
 */

input = getDirectory("Input directory");
output = getDirectory("Output directory");

Dialog.create("File type");
Dialog.addString("File suffix: ", ".TIF", 5);
Dialog.show();
suffix = Dialog.getString();

processFolder(input);

function processFolder(input) {
	list = getFileList(input);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + list[i]))
			processFolder("" + input + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
	// added by Ke
	list = getFileList(output);
		for (i = 0; i < list.length; i++) {
		    open(output + list[i]);
	}
	run("Concatenate...", "all_open title=[Concatenated Stacks]");
	saveAs("Tiff","H:\\GC Project\\RVKC314\\Concatenated Stacks.tif");
	
	
}

function processFile(input, output, file) {
	// do the processing here by replacing
	// the following two lines by your own code
	open(input + file);
	run("Grouped Z Project...", "projection=[Average Intensity] group=5");
	saveAs("Tiff", output + file);
	close();
	close();
	print("Processing: " + input + file);
	print("Saving to: " + output);
}
