# vim: ts=2 sw=2 noet

# The Rocker project makes container images for R.
# We'll be bringing in a specific Shiny version, so we can't use their Shiny
# container images.  Instead, we'll use their reproducable images.
# NOTE: This is where the R version is set!
FROM rocker/r-ver:4.4.2
LABEL org.opencontainers.image.base.name="rocker/r-ver:4.4.2"

# Set some image metadata!
# NOTE: Some of the OCI labels are set automatically by the GitHub Actions workflow.
#       For labes not set automatically, you'll need to set them both here
#       *and* in the workflow.
# NOTE: org.opencontainers.image.base.name is set above this block, as a subtle
#       reminder that, when you change the source image, you also need to
#       update the label!
LABEL org.opencontainers.image.documentation="TBD"
#LABEL org.opencontainers.image.description gets set to the repo's description
#LABEL org.opencontainers.image.licenses gets set to the license GitHub detects
#LABEL org.opencontainers.image.source gets set to the repo's URL
#LABEL org.opencontainers.image.title gets set to the repo's name
#LABEL org.opencontainers.image.url gets set to the repo's URL
#LABEL org.opencontainers.image.version gets set to the repo's branch/tag name

# Start with an apt-get upgrade, to update Ubuntu packages
RUN apt-get update && apt-get upgrade -y && apt-get clean && rm -rf /tmp/*

# Add the Microsoft ODBC Driver 18 package repo
# NOTE: If we change the upstream source to a different Ubuntu, we need to
# change this command!
RUN apt-get update && apt-get install -y curl && \
  curl -o /etc/apt/trusted.gpg.d/microsoft.asc \
  https://packages.microsoft.com/keys/microsoft.asc && \
  curl -o /etc/apt/sources.list.d/mssql-release.list \
  https://packages.microsoft.com/config/ubuntu/22.04/prod.list && \
	apt-get purge -y curl && apt-get autoremove -y

# Install any OS-level packages that we need for our R packages.
# unixodbc-dev is used by the odbc R package.
# zlib1g-dev is used by the httpuv R package.
# pkg-config is used since we're bringing in dev packages.
# ACCEPT_EULA is used by msodbcsql17
# NOTE: We avoid the driver for SQL Server 18, becaue it enables encryption by
# default.
ARG ACCEPT_EULA=y
RUN apt-get update && \
  apt-get install -y libodbc2 libcurl4-openssl-dev libssl-dev libxml2-dev msodbcsql17 pkg-config unixodbc-dev zlib1g-dev && \
	apt-get clean && rm -rf /tmp/*

# Copy all of the files into the container's scripts directory, and set that to
# be the working directory (the PWD) for all the R commands we run.
COPY . /usr/local/src/myscripts
WORKDIR /usr/local/src/myscripts

# Install renv, so that we can import the environment.
RUN R -q -e 'install.packages("renv")' && rm -rf /tmp/*

# Restore the environment from the lock file.
# NOTE: If you change the lock file, be sure to update the Ubuntu package list
# appropriately!
# NOTE: Renv caches downloaded packages, so we have to clean them up at the end.
#       The `use.cache` setting applies to installed packages, not sources.
# TODO: See if we can save the package downloads somewhere
RUN R -q -e 'renv::settings$use.cache(FALSE)' -e 'renv::restore()' && \
	rm -rf /tmp/* /root/.cache/R/renv/source/repository/*

# Declare the environment variables we want, and set some default values.
# NOTE: It's understood that these defaults won't work in production.
ENV SERVER=127.0.0.1
ENV DB=noDB
ENV DB_USERNAME=noUsername
ENV DB_PW=noPassword
ENV TZ=America/Los_Angeles
ENV DRIVER="ODBC Driver 17 for SQL Server"
ENV INSTITUTION_ID="SUNetID"

# When the container is run without an explicit command, this is what we do:
# Start our Shiny app!  Listen on port 3838, and expose that to the outside.
EXPOSE 3838
CMD ["R", "-q", "-e", "shiny::runApp('/usr/local/src/myscripts', host='0.0.0.0', port=3838)"]
