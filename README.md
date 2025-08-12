# Welcome to Arribere lab's github repository

### How to set up Git on Linux:
- Follow this [link](https://git-scm.com/downloads) to download the latest version of Git
- To set a Git username, open terminal and run the following command substituting with the username of your choice:
    `git config --global user.name "annie"`
- Double check this with this command, which should return the username you chose:
    `git config --global user.name` 
- For more information, follow this [link](https://help.github.com/en/github/using-git/getting-started-with-git-and-github)

### How to set up this repository on Linux:
- Open terminal and navigate to your home directory:
    `cd ~`
- Clone the arriberelab repository:
    `git clone https://github.com/arriberelab/arriberelab.git`
- Navigate into the cloned arriberelab directory and confirm it's there and contains the same files on GitHub:
    `cd arriberelab`
    `ls -la`
- For more information, follow this [link](https://help.github.com/en/github/getting-started-with-github/set-up-git#next-steps-authenticating-with-github-from-git) 

### How to save your GitHub password in Git:
- In terminal, enter the following:
    `git config --global credential.helper cache`
- Run the following command to adjust how long you want to store your password, where timeout is in seconds: 
    `git config --global credential.helper 'cache --timeout=3600' 
- Or if you want to store your Git credentials permanently, run:
    `git config credential.helper store`
- For more information, follow this [link](https://help.github.com/en/github/using-git/caching-your-github-password-in-git)

### How to push you changes to the arriberelab GitHub repository: 
- Ensure you're on the master branch by typing: 
    `git branch`
- If you are not on master, you can switch to it by running: 
    `git checkout master`
> Note: if you have uncommitted changes on your current branch, you must either [commit them](https://git-scm.com/docs/git-commit), [remove them](https://git-scm.com/docs/git-reset), or [stash them](https://git-scm.com/docs/git-stash) 
- Create/remove/edit a file in your local clone of the arriberelab repository: 
    `echo "print(\"hello world\")" > ~/arriberelab/hello_world.py`
- Track the modifciation in git with: 
    `git add ~/arriberelab/hello_world.py`
- Or if you've made several modifications, you can simply run: 
    `git add ~/arriberelab/*`
- Commit your changes: 
    `git commit -m "Add hello_world.py"`
- Push you changes to GitHub: 
    `git push`
- For more information on how git works follow this great introductory [tutorial](https://rogerdudler.github.io/git-guide/)
