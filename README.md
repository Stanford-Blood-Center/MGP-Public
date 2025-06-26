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
Note: MGP is currently only compatible with mTilda. MGP integration with 
other Laboratory Information Management Systems (LIMS) will require code changes.

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

MGP uses a number of APIs.  You do not need credentials to access these APIs,
but you do need to use a machine that can connect to the Internet.

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

  * On Linux and macOS, libodbc from the [unixODBC
    project](https://www.unixodbc.org).  Read more in the *ODBC* section, below.

  * The Microsoft [ODBC Driver for SQL
    Server](https://learn.microsoft.com/en-us/sql/connect/odbc/download-odbc-driver-for-sql-server?view=sql-server-ver17),
    **version 17**.

  * [zlib](https://zlib.net) (also known as "zlib1" or "zlib1g").

## ODBC

MGP uses ODBC to connect to mTilda's database.  ODBC separates its functionality
into the "driver manager"—which is database-agnostic—and a database-specific
driver.  In R, the [odbc
package](https://cran.r-project.org/web/packages/odbc/index.html) provides an
interface between the ODBC driver manager and R's [DBI
package](https://dbi.r-dbi.org).  The R packages (obbc and DBI) are installed
through renv; you are responsible for providing the driver manager and the
database driver.

On Windows, the ODBC driver manager is built-in.  On macOS and Linux, you will
have to install the [unixODBC project](https://www.unixodbc.org),

If you are using SQL Server to host the mTilda database, then your database
driver will be the Microsoft [ODBC Driver for SQL
Server](https://learn.microsoft.com/en-us/sql/connect/odbc/download-odbc-driver-for-sql-server?view=sql-server-ver17).
Note that **you must use Driver 17**, as Driver 18's defaults cause issues.

## Reference Databases

The `alignments.rda` file in the `ref` directories come from the ["Alignments"
folder](https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments) of the
[IPD-IMGT/HLA Database on GitHub](https://github.com/ANHIG/IMGTHLA/).  The
`updateAlignmentVersion.R` R script (located in the `functions` directory) is
responsible for building this file.  It is built in to the application to speed
up program load times, particularly for the Docker container (which cannot
cache files between executions).

The R script calls BLAASD (Build Loci Amino Acid Specific Dataframe), which
extracts alignment sequence information for a given locus from the
ANHIG/IMGTHLA GitHub repository to produce a dataframe of individual amino acid
positions for all alleles, for a user-defined HLA locus or loci.

The code block for BLAASD is derived from the HLAtools R package. Future 
versions of MGP may directly call HLAtools::buildAlignments().

We acknowledge the Mack Laboratory at UCSF for the usage of BLAASD:

Livia Tran, Ryan Nickens, Leamon Crooms IV, Derek Pappas, Vinh Luu, Josh Bredeweg,
Steven Mack, HLAtools: Toolkit for HLA Immunogenomics.
https://CRAN.R-project.org/package=HLAtools


We acknowledge the following publications describing the IPG-IMGT/HLA Database:

Dominic J Barker, Giuseppe Maccari, Xenia Georgiou, Michael A Cooper, Paul
Flicek, James Robinson, Steven G E Marsh, The IPD-IMGT/HLA Database, *Nucleic
Acids Research*, Volume 51, Issue D1, 6 January 2023, Pages D1053–D1060,
[https://doi.org/10.1093/nar/gkac1011](https://doi.org/10.1093/nar/gkac1011).

James Robinson, Dominic J Barker, Steven G E Marsh, 25 years of the
IPD-IMGT/HLA Database, *HLA*, Volume 103, Issue 6, June 2024, e15549,
[https://doi.org/10.1111/tan.15549](https://doi.org/10.1111/tan.15549).

J. Robinson, A. Malik, P. Parham, J. G. Bodmer, S. G. E. Marsh, IMGT/HLA - a
sequence database for the human major histocompatibility complex, *Tissue
Antigens*, March 2000, Volume 55, Issue 3, Pages 280-287,
[https://doi.org/10.1034/j.1399-0039.2000.550314.x](https://doi.org/10.1034/j.1399-0039.2000.550314.x).

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
The Dockerfile then installs the non-R libraries, installs renv, and hands off
to renv for R package installation.

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

Here is a demo of building the container image on a macOS system:

[![asciicast](https://asciinema.org/a/5TE4WiUwO1UtLulfU635s4BL1.svg)](https://asciinema.org/a/5TE4WiUwO1UtLulfU635s4BL1)

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

3. Download (or clone) a copy of this Git repository.  `cd` to your copy.

4. Run `renv restore`.  This will download and build the R pacakges.  Expect it
   to take some time.

You should now have an R environment that is ready to run MGP.

# How to Run

Once you have "built" MGP, you need to run it!  To do so, you first need to
identify which ODBC driver you will be using.  Then, you need to set
environment variables.  Finally, you can run MGP according to how you built it.

*For multi-user use, we recommend using
[ShinyProxy](https://www.shinyproxy.io)*: For MGP to work, it needs read access
to the mTilda database.  That means MGP needs access to database credentials,
and those should not be kept on end-user machines.  Running the MGP Docker
container inside ShinyProxy ensures the application code and configuration are
kept off of end-user machines.  In addition, ShinyProxy likely supports your
site's authentication system, and is also able to use "Social Auth" (like
Google) as a fallback.

## ODBC Driver

To configure MGP, you need to tell it which ODBC driver to use.
If you are using the Docker container, you can skip this step.

The ODBC driver manager provides a list of drivers, either in an app or in an
INI file.

For example, here is the contents of the `odbcinst.ini` file on a macOS system,
after installing the Microsoft ODBC Driver for SQL Server, version 17:

```
$ cat /opt/local/etc/odbcinst.ini
[ODBC Driver 17 for SQL Server]
Description=Microsoft ODBC Driver 17 for SQL Server
Driver=/opt/local/lib/libmsodbcsql.17.dylib
UsageCount=1
```

On Linux, the same file will typically be located at path `/etc/odbcinst.ini`.
On Windows, the information is available through the *ODBC Data Source
Administrator* application, in the *Drivers* tab.

![The "Drivers" tab of the Windows "ODBC Data Sourace Administrator"
application](https://github.com/user-attachments/assets/913a78ef-f008-432a-bb3f-b327b91fe513)

You will need to identify the name of the ODBC driver to use, and then set the
appropriate environment variable.

## Environment Variables

MGP is configured using environment variables.  The exact method of setting
environment variables depends on your operating system and how you are going to
run MGP.

Here are the environment variables you need to set:

* Name: `DRIVER`

  Derscription: The name of the ODBC Driver to use for connecting to the
  database.  This should be the name of a driver from the `odbcinst.ini`
  configuration file.  For example, "ODBC Driver 17 for SQL Server".

  If you are using the Docker container, this environment variable is already
  set for you in the container, so you should ignore it.

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

* Name: `MAINTAINER_EMAIL`

  Description: The email address of the person who maintains this installation
  of MGP.

  If you do not set the `MAINTAINER_EMAIL` environment variable, the string
  `maintaneremail@placeholder.com` will be used instead.

* Name: `INSTITUTION_ID`

  Description: The term that your site uses to refer to usernames.  At Stanford
  University, this is "SUNetID".

  If you do not set the `INSTITUTION_ID` environment variable, the string
  "Instituion ID" will be used instead.

The `TX`, `MAINTAINER_EMAIL`, and `INSTITUTION_ID` environment variables are
optional; all others are required.

## Running with Docker

If you used the *Docker* procedure to build MGP, you should have a container
image named "mgp:latest".

Before running the container, you will need to set environment variables.  You
can set these in your shell's local environment, in a text file of environment
variables, or on the `docker` command line.  We recommend placing environment
variables in a text file.

Here is an example file, named `mgp.env`, containing MGP environment variables:

```
SERVER=dbserver.example.com\mTilda
DB=mTilda
DB_USERNAME=something
DB_PW=somethingElse
MAINTAINER_EMAIL=admin@example.com
INSTITUTION_ID=MyBCID
TZ=America/New_York
```

(Since you are using the container, you should not set the `DRIVER` environment
variable.  The other variables should be customized to your environment.)

With the `mgp.env` file ready, you can use this Docker command to run the
container:

```
docker run --detach --rm --publish-all --env-file mgp.env --platform linux/amd64 mgp:latest
```

(Again, `--platform linux/amd64` is only needed if running Docker on a macOS
system.)

The MGP application listens on port 3838 inside the container; the
`--publish-all` operation will connect that port to a random port on your
machine.  The `docker container ls` command will tell you which port was
assigned.

The `--detach` option runs the container in the background, and the `--rm`
option will clean up the container after it is stopped.

As MGP runs, you can run the `docker logs` command to see the logs from the
container.  To stop the container, run `docker stop`.  For each Match Grade
that is performed, MGP will also generate a log file, which can be downloaded
through the MGP web site.

Here is an example of how to run the MGP Docker container:

[![asciicast](https://asciinema.org/a/724444.svg)](https://asciinema.org/a/724444)

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

Make sure you have set the environment variables.  Since you are not using the
Docker container, you *will* need to set the `DRIVER` environment variable.

To run the application, run the command `renv run app.R`.  After a short delay,
you should see the message `Listening on http://127.0.0.1:7346` (though you may
see a different number instead of "7346").  Follow the URL, and you should be
presented with the MGP main page.

As the program runs, it will print log messages to your terminal's standard
output.  For each Match Grade that is performed, MGP will also
generate a log file.

## Running in ShinyProxy

To use MGP within ShinyProxy, you should first confirm that the Docker
container works.  To do this, follow the instructions in the *Running with
Docker* section.  At the end, you should have a `mgp.env` file that configures
MGP for your environment.

Next, follow the ShinyProxy [Getting
Started](https://www.shinyproxy.io/documentation/getting-started/) instructions
to install ShinyProxy and perform basic configuration.  You may also want to
read the
[Configuration](https://www.shinyproxy.io/documentation/configuration/) section
for information on how to configure ShinyProxy to use your site's
authentication.  You should be able to use the pre-configured "Hello
Application" and "06\_tabsets" apps.

With ShinyProxy working, you may proceed to add MGP.  To do so, add a new app
to the list.  Here is an example configuration snippet:

```
proxy:
  docker:
    image-pull-policy: always
  specs:
    - id: mgp
      display-name: MGP
      description: Match Grade Populator
      container-image: ghcr.io/stanford-blood-center/mgp-public:main
      container-env-file: mgp.env
```

The above example uses the `mgp.env` file you created while following the
*Running with Docker* instructions.  It also sets ShinyProxy's image pull
policy to "always", which will prompt ShinyProxy to check for an updated
container image whenever an app is launched.

After changing the ShinyProxy configuration, restart ShinyProxy, and you should
have MGP available as an app:

![ShinyProxy's home page, showing the MGP
app](https://github.com/user-attachments/assets/2ca3f634-3dd6-43ce-9051-18c3cfd56a05)


# Copyright & Licensing

Match Grade Populator is © Stanford Blood Center, LLC.

The copyrightable parts of the IPD-IMGT/HLA database are covered under the
Creative Commons Attribution-NoDerivs License.  Read the [IPD-IMGT/HLA License
and Disclaimer](https://www.ebi.ac.uk/ipd/imgt/hla/licence/).
