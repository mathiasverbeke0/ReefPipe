# ReefPipe

ReefPipe is a command line workflow tool developed to address the challenging task of analyzing large amounts of metabarcoding data. It specifically focuses on processing COI and ITS metabarcoding data, automating parallel Amplicon Sequence Variant (ASV) inference and multi-reference taxonomic classification.

## Installation and Setup

To use ReefPipe, follow these setup steps. If you encounter installation issues, consider using the Docker image. It includes all required software components and can be accessed at [ReefPipe-Docker](https://github.com/mathiasverbeke0/ReefPipe-Docker). Using the Docker image eliminates individual installations and potential compatibility problems, as all dependencies are preconfigured.

### R Software

Make sure you have R version 4.1.3 or higher installed on your system. If you don't have it already, you can download and install it from the [official R project website](https://www.r-project.org/).

### Python Software

Make sure you have Python version 3.6 or higher installed on your system. If you don't have it already, you can download and install it from the [official Python website](https://www.python.org/).

### R Packages

R packages are automatically installed by ReefPipe scripts, eliminating the need for manual installation. However, errors might occur during the download of the DADA2 package, particularly on certain systems. In such cases, you may need to install additional software specific to your operating system. Please note that we cannot provide troubleshooting support for all operating systems, so it is recommended to address these issues independently.

### Python Modules

After installing Python, open your terminal and run the following command to install the necessary modules:

```bash
pip install boldigger_cline cutadapt argparse biopython tqdm pandas
```

### Clustal Omega Software

For Linux distributions that utilize the dnf package manager, such as Fedora, you can download and install Clustal Omega by executing the following command in your terminal: 

```bash
sudo dnf install clustal-omega
```

For Linux distributions that utilize the apt package manager, such as Ubuntu, you can download and install Clustal Omega by executing the following command in your terminal: 

```bash
sudo apt-get install clustalo
```

Please note that these examples cover 2 popular Linux distributions, but if you are using a different version or distribution, you will need to determine the appropriate command yourself. The package manager and installation method can vary across different Linux distributions.

For Windows and macOS users, a Clustal Omega binary has been incorporated into the ReefPipe script<sup>1</sup>.

### ReefPipe

If you have Git initialized on your system, you can clone or download the ReefPipe code from the GitHub repository by running the following command in your terminal:

```
git clone git@github.com:mathiasverbeke0/ReefPipe.git
```

If you don't have Git initialized on your system, you can download the ReefPipe code as a ZIP file by following these steps:

1. Visit the ReefPipe GitHub repository at https://github.com/mathiasverbeke0/ReefPipe.
2. Click on the green "Code" button.
3. Select "Download ZIP" from the dropdown menu.
4. Once the ZIP file is downloaded, extract its contents to a location of your choice on your system.

## Documentation
For detailed instructions on how to use ReefPipe, refer to the documentation.

## References
Please find below the list of references that are relevant to ReefPipe:

<sup>1</sup> Example Reference 1

Please refer to these references for further reading and to understand the underlying concepts and methodologies used in ReefPipe.