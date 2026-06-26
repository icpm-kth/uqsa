# The PATH variable

This is an explanation about the PATH variable and has a very distant
relation to the `uqsa` package.

It is just generally true about shells (e.g. `bash`) and can be useful
when troubleshooting the model-building procedures. Consult this text
when a *command is not found*.

## Locations of Applications

Most unix like systems (\*BSD, GNU/Linux, macOS) have a way to find
programs: the environment variable called `PATH`; a colon (`:`)
separated string.

This is a possible value:

``` sh
echo "$PATH"
#> /home/andreikr/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/snap/bin
```

If you need to actually read it:

``` sh
echo "$PATH" | tr ':' '\n'
#> /home/andreikr/.local/bin
#> /usr/local/sbin
#> /usr/local/bin
#> /usr/sbin
#> /usr/bin
#> /sbin
#> /bin
#> /usr/games
#> /usr/local/games
#> /snap/bin
#> /snap/bin
```

The package manager of unix-like systems installs applications into
places that are listed in the `PATH` variable. This way they can be
found by name, rather than typing out the entire path.

``` sh
(
cat<<EOF
#!/bin/sh
echo "Hello World"
EOF
) > test.sh
chmod 744 test.sh
[ -f test.sh ] && test.sh
./test.sh && rm test.sh
#> sh: 8: test.sh: not found
#> Hello World
```

This code-block creates a file (`test.sh`), makes it executable, and
then tries to execute it. This doesn’t work, even though we are in the
same directory as the file. But, typing out a relative path `./test.sh`,
using `.`, it *is* found.

### Built-in commands

Some commands are shell *built-ins* (e.g.: `echo`, `printf`, `pwd`),
they don’t need the `PATH` variable at all as they are an integral part
of the shell itself (different shells can have slightly different
built-ins).

Other commands *are* distinct, executable files that need to be found
(e.g.: `env`, `awk`, `perl`). They can be symbolic links to such files
as well.

Here is an example for `awk`:

``` sh
type awk
which awk | xargs ls -gG
readlink -f `which awk`
#> awk is /usr/bin/awk
#> lrwxrwxrwx 1 21 Apr  8  2024 /usr/bin/awk -> /etc/alternatives/awk
#> /usr/bin/gawk
```

So, `awk` is not a built-in, but a file. The command `which` can find
it, so it is in `PATH`. It is also a *symbolic link*, we can follow that
link using `readlink -f` and find out that this particular system is
using *GNU awk* (`gawk`).

When a program is not in one of the directories in the PATH variable, it
can only be called by typing the entire path to it, or using a relative
path.

## Add locations to PATH

An interactive shell has a startup file:

- `~/.bashrc` for bash
- `~/.zshrc` for zsh (among others)

They also read `~/.profile`, which often has something like this in it
(if it doesn’t, you can add this yourself):

``` sh
# set PATH so it includes user's private bin if it exists
if [ -d "$HOME/.local/bin" ] ; then
    PATH="$HOME/.local/bin:$PATH"
fi
```

So, if `~/.local/bin` exists, it is appended to `PATH`, and you can
definitely use that location to store shell scripts (just create it).

``` sh
mkdir -p ~/.local/bin
```

There is a lot of fine detail related to the concepts of *login shell*
or *(non-) interactive shell* and what shell reads which files in either
of those cases (this is a bit complicated).

### IFF still not in PATH

A new entry can be added from within the `~/*rc` files, e.g.:

``` sh
## for bash
[ "$SHELL" == "/bin/bash" ] && echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
## for zsh
[ "$SHELL" == "/bin/zsh" ] && echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc
```
