# Match Grade Populator

Match Grade Populator (MGP) is a novel computational tool that automates the 
calculations used in bone marrow transplant compatibility evaluation.

MGP is written in R, and uses the Shiny web application framework.

# What it does

MGP calculates the number of HLA phenotype mismatches, quantifies the number of 
Graft vs Host (GvH) and Host vs Graft (HvG) directional mismatches, evaluates 
donor compatibility against allogeneic HLA antibodies, and assesses T-Cell 
Epitope permissibility for DPB1 mismatched alleles.

# How it works

On start up, MGP prompts the user for a mTilda username upon start-up. 

1. Enter a registered mTilda username and press '**Validate username**'. If the entered
name is not found, the user will not be able to access the application.
2. Enter a patient ITL in the text box labeled '**Patient ITL**'.
3. A drop-down box populated with donors available for evaluation for the entered
patient ITL. Select the desired donor for evaluation. Only one donor may be 
evaluated at a time.
4. Click on the '**Run MGP**' button to begin evaluation. A red, busy spinner will 
appear in the main panel while the program is running.

Output to the main panel is dynamically adjusted based on the mTilda user’s authorization level. 

Technologist level credentials trigger ‘**Calculation**’ mode, which writes calculated 
results to the Match Grade software. 

Higher tier credentials trigger ‘**Review**’ mode, which shows any  discrepancies 
between MGP calculations and Match Grade content. 

Both modes display locus-level mismatches, if donor serum antibodies are present, 
and whether surrogates were used due to donor alleles not being present in 
Class I and Class II Single Antigen Bead (SAB) panels. 

Logs containing important information on the evaluation process are available 
for download after each run. 


# APIs

#### **IPD-IMGT/HLA's DPB1 T-Cell Epitope Algorithm v 2.0**
- The HLA-DPB1 T-cell Epitope API determines HLA-DPB1 T-Cell epitope permissibility.
- Documentation for the API can be found [here](https://www.ebi.ac.uk/ipd/imgt/hla/matching/match_apis/).


#### **National Marrow Donor Program (NMDP) Multiple Allele Codes (MAC)**
- The NMDP MAC API returns possible subtypes for a given NMDP allele.
- Documentation for the API can be found [here](https://hml.nmdp.org/mac/).
- **Note**: if the API is down, MGP will fail for any evaluations that consist of NMDP
genotypes. 

# Prerequisites

MGP is built on the following software:

* **R 4.4**: The latest version in the 4.4.x series should work.  

* **renv**: The [renv](https://rstudio.github.io/renv/) R package is used to
  download and install the other packages that MGP needs.  This is the only R
  package that you should need to download & install manually.

The R packages we use rely on the following libraries:

  * [OpenSSL](https://openssl-library.org) (particularly `libssl`)

  * Curl (particularly [libcurl](https://curl.se/libcurl/), configured for
    OpenSSL)

  * [libxml2](https://gitlab.gnome.org/GNOME/libxml2/-/wikis/home)

  * libodbc from the [unixODBC project](https://www.unixodbc.org)

  * The Microsoft [ODBC Driver for SQL
    Server](https://learn.microsoft.com/en-us/sql/connect/odbc/download-odbc-driver-for-sql-server?view=sql-server-ver17),
    **version 17**.

  * [zlib](https://zlib.net) (also known as "zlib1" or "zlib1g").

## Reference Databases

In addition to software packages, MGP relies on several data sources.

TODO: Add text covering the remaining items in the `ref` directory.

# How to Build

Even though MGP is an R app, and does not need to be "compiled" in the way that
a Java program is compiled, the app still must be "built": It has a number of
dependencies—Shiny, for example—which must be downloaded and installed.  At
least one of those dependencies (the `odbc` package) depends on an ODBC library.
And R modules themselves are often compiled.

There are three ways to build MGP.  From easiest to hardest the options are…

1. You can use Docker

2. You can use GitHub Actions

3. You can build from source.

## Docker

The easiest way of building MGP is to use the Dockerfile.  The container
image uses base images from the [Rocker Project](https://rocker-project.org),
which provides a R installation built on top of Ubuntu.
The Dockerfile then installs the libraries, renv, and hands off to renv for R
package installation.

To build the docker image, follow these steps:

1. Download & install Docker (or a compatible Docker-alike).

2. Download (or clone) a copy of this Git repository.

3. Run the following command: `docker buildx build --platform linux/amd64 --tag mgp .`

   The `--platform linux/amd64` portion of the command ensures that the build
   runs on 64-bit Linux on Intel/AMD CPUs, even if you are running the build
   from a Mac.  Native builds on the ARM platform are not available right now.

After the container build completes, you will have a container image named
(or "tagged") "mgp".  More specifically, it will be tagged "mgp:latest", as
this is your latest build of the MGP container image.

## GitHub Actions

This repository includes a GitHub Actions workflow which will run any time
there is a push.  The workflow builds the container image, making it available
in the repository's *Packages* area.  Click on the package name to view
instructions on how to pull the container image.

You should expect your container image's name to be based on the GitHub repo
and branch names.  For example, if your repo is named "matchgrade" and the
branch name is "sbc", the container image will be tagged "matchgrade:sbc".

## Building from Source

Building MGP from source is only really useful if you plan on doing development
work on MGP.  If you are just trying to use MGP, then building from source is
probably not worth it.

These instructions will be fairly generic.  Converting the generic instructions
to specific instructions is left as an exercise for the reader.

1. Download and install R and renv.

2. Install all of the libraries listed in the *Requirements* section.  Ensure
   that both the library *and the header files* are available.

3. TODO: ODBC Configuration

4. Download (or clone) a copy of this Git repository.  `cd` to your copy.

4. Run `renv restore`.  This will download and build the R pacakges.  Expect it
   to take some time.

You should now have an R environment that is ready to run MGP.

# How to Run

Once you have "built" MGP, you need to run it!  To do so, you first need to set
environment variables.  Then, you need to run MGP according to how you built
it.

## Environment Variables

MGP is configured using environment variables.  The exact method of setting
environment variables depends on your operating system and how you are going to
run MGP.

Here are the environment variables you need to set:

* Name: `SERVER`

  Description: The hostname or IP address of the SQL Server.

* Name: `DB`

  Description: The name of the SQL Server database.

* Name: `DB_USERNAME`

  Description: A username with read-only access to the database.

* Name: `DB_PW`

  Description: The password for `DB_USERNAME`.

* Name: `TZ`

  Description: The time zone to use for log messages.

  For the `TZ` environment variable, you should set it to a Canonical TZ
  identifier from the [List of tz database time
  zones](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones).  For
  example, if your machine uses Japan Standard Time, set the `TZ` environment
  variable to `Asia/Tokyo`.  If you try setting it to `Tokyo`, which is not a
  Canonical identifier, things will break.

  If you do not set the `TZ` environment variable, the behavior depends on how
  you run MGP: MGP may automatically use the system time zone, or it may fall
  back to UTC.
  
The `TZ` environment variable is the only one that is optional: All
other environment variables are required.

## Running with Docker

If you used the *Docker* procedure to build MGP, you should have a container
image named "mgp:latest".

TODO: ODBC Configuration.  Running.

As the program runs, it will print log messages to your terminal's standard
output.  Those messages can also be viewed using the `docker container logs`
command.  For each Match Grade that is performed, MGP will also
generate a log file, which can be downloaded through the MGP web site.

## Running with GitHub Actions

If you used the *GitHub Actions* procedure to build MGP, your repository should
contain a package, which is your container image.  The package name and image
tag will depend on your GitHub repository and branch names.  If you click on
the package name in GitHub, you will be presented with instructions on how to
download (or "pull") the package.  Follow those instructions.

After pulling the package, follow the instructions in *How to Run: Docker*.

## Running from Source

If you used the *Building from Source* procedure to build MGP, your Git
repository's `renv` directory should now contain lots of files, representing
the R packages that renv installed.

Make sure you have set the environment variables.

To run the application, run the command `renv run app.R`.  After a short delay,
you should see the message `Listening on http://127.0.0.1:7346` (though you may
see a different number instead of "7346").  Follow the URL, and you should be
presented with the MGP main page.

As the program runs, it will print log messages to your terminal's standard
output.  For each Match Grade that is performed, MGP will also
generate a log file.

# Copyright & Licensing

Match Grade Populator is © Stanford Blood Center, LLC.

TODO: Add copyright/licensing information about the items in the `ref`
directory.
