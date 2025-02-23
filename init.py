def initialize():
    import sys
    import subprocess

    # Install requirements from requirements.txt
    def install_requirements():
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])
            print("Requirements installed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Failed to install requirements: {e}")

    # Clone a repository from GitHub
    def clone_repo(repo_url, clone_dir):
        try:
            subprocess.check_call(["git", "clone", repo_url, clone_dir])
            print(f"Repository cloned successfully into {clone_dir}.")
        except subprocess.CalledProcessError as e:
            print(f"Failed to clone repository: {e}")

    install_requirements()
    clone_repo("https://github.com/vaskiax/celestials")