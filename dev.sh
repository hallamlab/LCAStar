HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
NAME=lcastar
UTILS=$HERE/src/$NAME/utils.py
python --version || exit 1
USER=$(python $UTILS USER)
VER=$(python $UTILS VERSION)
DOCKER_IMAGE=quay.io/$USER/$NAME

CONDA=conda
# CONDA=mamba # https://mamba.readthedocs.io/en/latest/mamba-installation.html#mamba-install
echo "$NAME:$VER"
echo ""

# this file contains a list of commands useful for dev,
# providing automation for some build tasks

# example workflow 1, pip:
# dev.sh --idev # create a local conda dev env
# # add pypi api token as file to ./secrets [https://pypi.org/help/#apitoken]
# # make some changes to source
# # bump up ./src/fabfos/version.txt
# dev.sh -bp # build the pip package
# dev.sh -up # test upload to testpypi
# dev.sh -upload-pypi # release to pypi index for pip install

# example workflow 2, conda:
# dev.sh --idev # create a local conda dev env
# dev.sh -bp # build the pip package
# dev.sh -bc # build conda package from pip package
# dev.sh -uc # publish to conda index

# example workflow 3, containerization:
# dev.sh --idev # create a local conda dev env
# dev.sh -bd # build docker image
# dev.sh -ud # publish to quay.io
# dev.sh -bs # build apptainer image from local docker image

case $1 in
    ###################################################
    # environments

    --idev) # with dev tools for packaging
        env_name=${NAME}_dev
        cd $HERE/envs
        echo "creating new conda env: $env_name"
        echo "WARNING: you will need to install docker and apptainer individually"
        sleep 2
        $CONDA env create --no-default-packages -n $env_name -f ./base.yml \
        && $CONDA env update -n $env_name -f ./dev.yml
    ;;
    --ibase) # base only
        env_name=${NAME}
        cd $HERE/envs
        echo "creating new conda env: $env_name"
        sleep 2
        $CONDA env create --no-default-packages -n $env_name -f ./base.yml
    ;;

    ###################################################
    # build

    -bp) # pip
        # build pip package
        rm -r build
        rm -r dist
        python -m build
    ;;
    -bpi) # pip - test install
        pip instalconda l $HERE/dist/$NAME-$VER-py3-none-any.whl
    ;;
    -bpx) # pip - remove package
        pip uninstall -y $NAME
    ;;
    -bc) # conda
        # requires built pip package
        rm -r $HERE/conda_build
        python ./conda_recipe/compile_recipe.py
        $HERE/conda_recipe/call_build.sh
    ;;
    -bd) # docker
        echo "$DOCKER_IMAGE"
        docker build -t $DOCKER_IMAGE:$VER .
    ;;
    -bs) # apptainer image *from docker*
        echo "$DOCKER_IMAGE"
        apptainer build $NAME.sif docker-daemon://$DOCKER_IMAGE:$VER
    ;;

    ###################################################
    # upload

    -up) # pip (testpypi)
        PYPI=testpypi
        TOKEN=$(cat secrets/${PYPI}) # https://pypi.org/help/#apitoken
        python -m twine upload --repository $PYPI dist/*.whl -u __token__ -p $TOKEN
    ;;
    -upload-pypi) # pip (pypi)
        echo "not all dependencies are available on pypi, so this is not a good idea..."
        # PYPI=pypi
        # TOKEN=$(cat secrets/${PYPI}) # https://pypi.org/help/#apitoken
        # python -m twine upload --repository $PYPI dist/*.whl -u __token__ -p $TOKEN
    ;;
    -uc) # conda (personal channel)
        # run `anaconda login` first
        find ./conda_build -name *.tar.bz2 | xargs -I % anaconda upload -u $DEV_USER %
    ;;
    -ud) # docker
        echo "$DOCKER_IMAGE"
        # login and push image to quay.io
        # sudo docker login quay.io
	    docker push $DOCKER_IMAGE:$VER
        echo "!!!"
        echo "remember to update the \"latest\" tag"
        echo "https://$DOCKER_IMAGE?tab=tags"
    ;;
    
    ###################################################
    # run

    -r)
        shift
        export PYTHONPATH=$HERE/src:$PYTHONPATH
        python -m $NAME $@
    ;;
    -rd) # docker
            # -e XDG_CACHE_HOME="/ws"\
        echo "$DOCKER_IMAGE"
        shift
        docker run -it --rm \
            -u $(id -u):$(id -g) \
            --mount type=bind,source="$HERE",target="/ws"\
            --workdir="/ws" \
            $DOCKER_IMAGE:$VER /bin/bash
    ;;
    -rs) # apptainer
            # -e XDG_CACHE_HOME="/ws"\
        shift
        apptainer exec \
            --bind ./:/ws \
            --workdir /ws \
            $HERE/$NAME.sif fabfos /bin/bash
    ;;

    ###################################################
    # docs

    -d)
        sphinx-build -M html ./docs/src ./docs/build
    ;;
    
    ###################################################
    # test

    -t) # single manual test
        echo "no test"
    ;;

    ###################################################

    *)
        echo "bad option"
        echo $1
    ;;
esac
