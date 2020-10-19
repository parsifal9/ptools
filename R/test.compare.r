test.compare <- function(try, truth, test.name, testing)
{
	truth.file <- paste("truth", "test", test.name, "r", sep = ".")

	if (!testing)
	{
		if (file.exists(truth.file))
			msg <- paste("Overwriting 'truth file' for test \"", test.name,
				"\".", sep = "")
		else
			msg <- paste("Writing new 'truth file' for test \"", test.name,
				"\".", sep = "")

		warning(msg)
		##dput(try, file=truth.file)
		truth <- try
		save("truth", file = truth.file)
		return(TRUE)
	}

	if (!file.exists(truth.file))
	{
		msg <- paste("Test failed because 'truth file', \"", truth.file,
			"\", is missing.", sep = "")
		warning(msg)
		return(FALSE)
	}

	##truth <- dget(truth.file)
	load(truth.file)
	diffs <- all.equal(try, truth)
	ok <- identical(diffs, TRUE)
	if (!ok)
	{
		msg <- paste("Failed test of \"", test.name, "\".", sep = "")
		warning(msg)
		return(list(try = try, truth = truth, diffs = diffs))
	}

	TRUE
}
