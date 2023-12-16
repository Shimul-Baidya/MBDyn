# Contributing Guidelines
There are many ways to contribute to MBDyn:

 - using it and reporting bugs, fixing/adding tutorials, documentation, and so on;
 - submitting patches that fix bugs or implement new features, and so on;
 - establishing a grant with DAER/Polimi to implement the features or develop 
   the models you need.

# Access to the Gitlab repository
The MBDyn main repository was recently moved to this Gitlab repository, 
hosted on 
[POLIMI](https://www.polimi.it/) servers. At the moment, access to some
Gitlab functionalities for external users is limited. 
~~we're working with 
the ITC department of [POLIMI](https://www.polimi.it/) to improve it. 
However, in the meantime,
**even for registered POLIMI or Google users it is necessary to be registered as Project Members**
to be able to fork the project and therefore also to submit merge requests.~~


Be advised, however, that we are not willing to grant easy project membership,
because it may have critical drawbacks; for example,
pushing random binary blobs into the Gitlab repository
(yes, this happened: someone committed executables and object files!) 
would slow everyone's work, and would require substantial cleanup effort on our side.

Thus, we are only willing to give you write access to the repository after we are sure that we can
trust your work. In the meantime, you can always work on your own fork and host it on any public GIT,
e.g., Github, Gitlab, or wherever it's best for you. And, of course, you can open
issues, ask for feedback and request us to merge your work into our develop branch.

While we work out a better solution, to be granted
access to the Gitlab repository please contact us through the 
[MBDyn users mailing list](https://www.mbdyn.org/Mailing.html).

# MBDyn Developers Guidelines
MBDyn is developed primarily by _internal_ developers at Politecnico di 
Milano - Department of Aerospace Science and Technology 
([DAER](http://www.aero.polimi.it/)).  
Contributions from _external_ developers are welcome, of course.
As the project's maintainers, we only ask interested coders to follow the 
simple steps here presented.
From the moment the main repository was moved to Git, we switched to the branching 
model described [here](https://nvie.com/posts/a-successful-git-branching-model/), 
so if you are a new developer please read that page carefully prior to committing
and pushing changes.

## Guidelines for _external_ developers
As a general rule, before starting to write your piece of code, consider opening an 
[issue](https://gitlab.com/help/user/project/issues/index.md), to discuss your intentions
with other developers. 
~~This, at the moment, **requires you to be
recognised as Project Member**, so please ask us to activate your account by 
writing to the [MBDyn users mailing list](https://www.mbdyn.org/?Mailing_Lists).~~

To open an issue, you only have to be signed in to the website: you can do it either with 
your Google or [POLIMI](https://www.polimi.it/) account.


Developers that would like to contribute to MBDyn but are not in the project
regulars must:
 - fork the Gitlab repository hosted [here](https://gitlab.polimi.it/Pub/mbdyn.git)
 - checkout a fresh branch from the `develop` branch
 - commit and push to your branch in your forked repository
 - issue a merge request into the develop branch of MBDyn main repository

## Guidelines for _internal_ developers
Developers that would like to be included among the regulars in MBDyn must:

 - request "developer" access permission to the MBDyn administrators, 
      @andomasarati in primis
 - clone the repository from [here](https://gitlab.polimi.it/Pub/mbdyn.git)
 - checkout a fresh branch from the `develop` branch
 - before pushing, make sure that your contribution follows the code development model 
      indicated above
 - whenever you feel that your contribution would benefit from a discussion
      with other developers, issue a merge request instead of directly pushing
      to the main repository

