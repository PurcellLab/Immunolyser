# Immunolyser 2.0 Local Version Setup Guide

**Why Docker?**\
Using Docker ensures that all the tools and dependencies are correctly set up in a container, so you don't need to worry about installing them individually on your system. It guarantees compatibility and simplifies the process of setting up Immunolyser 2.0.

### Notes:

- **Supported Systems**:\
  This guide supports the following operating systems:

  - **Windows** (10/11)
  - **macOS**
  - **Linux**

## Prerequisites

1. **Docker Desktop**\
   You’ll need Docker Desktop installed on your machine. Docker is used to run applications in containers, ensuring compatibility across systems.
   - **Download Docker Desktop**: [Docker Desktop Download](https://www.docker.com/products/docker-desktop) 
   - **System Requirements**: Windows (10/11), macOS, or Linux. You will need administrator privileges to run Docker.

2. **Disk Space**\
   Docker requires disk space for storing container images, volumes, and other related files. Make sure you have enough space available on your system:
   - For the **Immunolyser** Docker setup, the combined disk space used by the required images is approximately **29 GB**.
   - It’s recommended to have at least **50 GB** of free space to ensure smooth operation, as Docker images, container logs, and temporary files can take up additional space.
   - You can monitor the disk usage in Docker Desktop under the **Settings > Resources > Disk** tab.

## Setup Instructions

### 1. **Download Required Files**

You need to download the following files:

- **Dockerfile**: <a href="https://github.com/prmunday/Immunolyser/raw/develop/Dockerfile" target="_blank">Download Dockerfile</a>  
- **docker-compose.yml**: <a href="https://github.com/prmunday/Immunolyser/raw/develop/docker-compose.yml" target="_blank">Download docker-compose.yml</a>  
- **entrypoint.sh**: <a href="https://github.com/prmunday/Immunolyser/raw/develop/entrypoint.sh" target="_blank">Download entrypoint.sh</a>  

#### Saving the Files:
To save the files, right-click on each link and select **"Save Link As"**, or press **Ctrl+S** (Windows) / **Cmd+S** (Mac). Ensure that:
- `Dockerfile` is saved without a `.txt` extension.  
- `entrypoint.sh` is saved with the `.sh` extension.  
- `docker-compose.yml` is saved with the `.yml` extension.  

### 2. **Create Required Folders**

In the same directory where you downloaded the Dockerfile, create the following folders:

- `data_files`
- `data_files_2`
- `tools`

### 3. **Download and Place Tools**

To make the Docker setup work, you'll need to download several tools and place them in the `tools` folder. Follow the steps below:

- **Seq2Logo**:

  - <a href="https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=seq2logo&version=2.1&packageversion=2.1&platform=all" target="_blank">Download Seq2Logo 2.1</a>.
  - Download the `.gz` file and save it in the `tools` folder, next to the Dockerfile.

- **GibbsCluster 2.0f**:

  - <a href="https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=gibbscluster&version=2.0&packageversion=2.0f&platform=Linux" target="_blank">Download GibbsCluster 2.0f</a>.
  - Save the `.gz` file in the `tools` folder.

- **NetMHCpan 4.1b**:

  - <a href="https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCpan&version=4.2&packageversion=4.2b&platform=Linux" target="_blank">Download NetMHCpan 4.2b</a>.
  - Save the `.gz` file in the `tools` folder.
  - Ensure the file is named `netMHCpan-4.1b.Linux.tar.gz` (remove any extra numbers added by the system).

- **NetMHCIIpan 4.3f**:

  - <a href="https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=netMHCIIpan&version=4.3&packageversion=4.3i&platform=Linux" target="_blank">Download NetMHCIIpan 4.3i</a>.
  - Save the `.gz` file in the `tools` folder as `netMHCIIpan-4.3f.Linux.tar.gz`.

### **Folder Structure**  

After downloading the tools, your `tools` directory should be structured as shown below:

<p align="center">
  <img src="https://raw.githubusercontent.com/prmunday/Immunolyser/main/Documentation/Screenshots/Local%20Version%20Directory%20Screenshot.png" alt="Tools Folder Structure" />
</p>

### 4. **Build and Run Docker**

After you've set up the tools, you can build and run the Docker container. Follow these steps:

1. Open **Terminal** (macOS/Linux) or **Command Prompt** (Windows).
2. Navigate to the directory where the Dockerfile is located.
3. Run the following commands:
   ```sh
   docker-compose build
   ```

- Docker will take **10 to 30 minutes** to build the container, depending on your system.
- Make sure you have enough disk space for the build process.
4. Once the build is complete, run the following command to start the container:
   ```sh
   docker-compose up
   ```
5. Once the container is running, you can access the app at http://localhost:5001/.

### 5. **Check for Errors**

If any errors occur during the Docker setup, check that all the files are in place and named correctly.

### 6. **Managing Disk Space**

Each time you build a new Docker image, old images remain on your system, taking up disk space. To free up space, remove unused images and containers periodically by running:

```sh
docker system prune -a
```

This command will delete all stopped containers, networks, and dangling images. Use it cautiously as it removes all unused Docker objects.

---

- **Troubleshooting**:\
  If you encounter issues, check the [Docker Troubleshooting Guide](https://docs.docker.com/get-docker/) for help with installation or configuration.



