/*
 * David Galanter and John Collins
 */

/*
 * Global variables
 */
var meshResolution;

// Particle states
var mass;
var vertexPosition, vertexNormal, newPos, newVel;
var vertexVelocity;

// Spring properties
var K, restLength; 

// Force parameters
var Cd;
var uf, Cv;


/*
 * Getters and setters
 */
function getPosition(i, j) {
    var id = i*meshResolution + j;
    return vec3.create([vertexPosition[3*id], vertexPosition[3*id + 1], vertexPosition[3*id + 2]]);
}

function setPosition(i, j, x) {
    var id = i*meshResolution + j;
    vertexPosition[3*id] = x[0]; vertexPosition[3*id + 1] = x[1]; vertexPosition[3*id + 2] = x[2];
}

function getNormal(i, j) {
    var id = i*meshResolution + j;
    return vec3.create([vertexNormal[3*id], vertexNormal[3*id + 1], vertexNormal[3*id + 2]]);
}

function getVelocity(i, j) {
    var id = i*meshResolution + j;
    return vec3.create(vertexVelocity[id]);
}

function setVelocity(i, j, v) {
    var id = i*meshResolution + j;
    vertexVelocity[id] = vec3.create(v);
}


/*
 * Provided global functions (you do NOT have to modify them)
 */
function computeNormals() {
    var dx = [1, 1, 0, -1, -1, 0], dy = [0, 1, 1, 0, -1, -1];
    var e1, e2;
    var i, j, k = 0, t;
    for ( i = 0; i < meshResolution; ++i )
        for ( j = 0; j < meshResolution; ++j ) {
            var p0 = getPosition(i, j), norms = [];
            for ( t = 0; t < 6; ++t ) {
                var i1 = i + dy[t], j1 = j + dx[t];
                var i2 = i + dy[(t + 1) % 6], j2 = j + dx[(t + 1) % 6];
                if ( i1 >= 0 && i1 < meshResolution && j1 >= 0 && j1 < meshResolution &&
                     i2 >= 0 && i2 < meshResolution && j2 >= 0 && j2 < meshResolution ) {
                    e1 = vec3.subtract(getPosition(i1, j1), p0);
                    e2 = vec3.subtract(getPosition(i2, j2), p0);
                    norms.push(vec3.normalize(vec3.cross(e1, e2)));
                }
            }
            e1 = vec3.create();
            for ( t = 0; t < norms.length; ++t ) vec3.add(e1, norms[t]);
            vec3.normalize(e1);
            vertexNormal[3*k] = e1[0];
            vertexNormal[3*k + 1] = e1[1];
            vertexNormal[3*k + 2] = e1[2];
            ++k;
        }
}

var clothIndex, clothWireIndex;
function initMesh() {
    var i, j, k;

    vertexPosition = new Array(meshResolution*meshResolution*3);
    vertexNormal = new Array(meshResolution*meshResolution*3);
    newPos = new Array(meshResolution*meshResolution);
    newVel = new Array(meshResolution*meshResolution);

    clothIndex = new Array((meshResolution - 1)*(meshResolution - 1)*6);
    clothWireIndex = [];

    vertexVelocity = new Array(meshResolution*meshResolution);
    restLength[0] = 4.0/(meshResolution - 1);
    restLength[1] = Math.sqrt(2.0)*4.0/(meshResolution - 1);
    restLength[2] = 2.0*restLength[0];

    for ( i = 0; i < meshResolution; ++i )
        for ( j = 0; j < meshResolution; ++j ) {
            setPosition(i, j, [-2.0 + 4.0*j/(meshResolution - 1), -2.0 + 4.0*i/(meshResolution - 1), 0.0]);
            setVelocity(i, j, vec3.create());

            if ( j < meshResolution - 1 )
                clothWireIndex.push(i*meshResolution + j, i*meshResolution + j + 1);
            if ( i < meshResolution - 1 )
                clothWireIndex.push(i*meshResolution + j, (i + 1)*meshResolution + j);
            if ( i < meshResolution - 1 && j < meshResolution - 1 )
                clothWireIndex.push(i*meshResolution + j, (i + 1)*meshResolution + j + 1);
        }
    computeNormals();

    k = 0;
    for ( i = 0; i < meshResolution - 1; ++i )
        for ( j = 0; j < meshResolution - 1; ++j ) {
            clothIndex[6*k] = i*meshResolution + j;
            clothIndex[6*k + 1] = i*meshResolution + j + 1;
            clothIndex[6*k + 2] = (i + 1)*meshResolution + j + 1;
            clothIndex[6*k + 3] = i*meshResolution + j;
            clothIndex[6*k + 4] = (i + 1)*meshResolution + j + 1;            
            clothIndex[6*k + 5] = (i + 1)*meshResolution + j;
            ++k;
        }
}


/*
 * KEY function: simulate one time-step using Euler's method
 */


// s is 0,1,2 represents type of spring
function getSpringForce(p, q, s){
    var diff = vec3.subtract(vec3.create(p),q);
    var len = vec3.length(diff);    
    return vec3.scale(diff, K[s]*(restLength[s]-len)/len);
    //return vec3.scale(diff, 1000.0*(restLength[s]-len)/len);
}

function simulate(stepSize) {
    for(i = 0; i < meshResolution*meshResolution; ++i){
        newPos[i] = vec3.create();
        newVel[i] = vec3.create();
    }

    var end = meshResolution - 1;
    for( y = 0; y < meshResolution; ++y){
        for( x = 0; x < meshResolution; ++x){
            var force = vec3.create();  // force for this particle in mesh

            // structural
            if(x > 0){  // left
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y,x-1),0));
            }
            if(x < end){ // right
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y,x+1),0));
            }
            if(y > 0){  // down
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y-1,x),0));
            }
            if(y < end){ // up
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y+1,x),0));
            }
            
            // shear
            if(x > 0 && y > 0){ // down left
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y-1,x-1),1));
            }
            if(x < end && y > 0){ // down right
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y-1,x+1),1));
            }
            if(x > 0 && y < end){ // up left
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y+1,x-1),1));
            }
            if(x < end && y < end){ // up right
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y+1,x+1),1));
            }
            
            // flexion
            if(x > 1){  // left
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y,x-2),2));
            }
            if(x < end - 1){ // right
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y,x+2),2));
            }
            if(y > 1){  // down
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y-2,x),2));
            }
            if(y < end - 1){ // up
                vec3.add(force,getSpringForce(getPosition(y,x),getPosition(y+2,x),2));
            }
            
            // gravity
            force[1] = force[1] - mass * 9.8;
            
            // damping
            vec3.add(force, vec3.scale(getVelocity(y,x),-Cd));
            
            // viscous
            var fc = vec3.dot(getNormal(y,x), vec3.subtract(vec3.create(uf),getVelocity(y,x)));
            vec3.add(force, vec3.scale(getNormal(y,x), Cv * fc));
            
            // calculate and save new position and velocity for this particle
            var id = y * meshResolution + x;
            if(y == end && (x == 0 || x == end)){    // if at corners just set new to current position
                newPos[id] = getPosition(y,x);
            }else{
                //HERE'S THE BEGINNING OF THE CHANGE JOHNNAAAAAY
                //newPos[id] = vec3.add(getPosition(y,x),vec3.scale(getVelocity(y,x),stepSize));
                newVel[id] = vec3.add(getVelocity(y,x),vec3.scale(force, stepSize/mass));
                newPos[id] = vec3.add(getPosition(y,x),vec3.scale(vec3.create(newVel[id]),stepSize));
                //HERE'S THE END OF THE CHANGE JOHNNAAAAAY
            }

        }
    }

    // now actually set the values
    for(y = 0; y < meshResolution; ++y){
        for(x = 0; x < meshResolution; ++x){
            var id = y * meshResolution + x;
            setPosition(y,x,newPos[id]);
            setVelocity(y,x,newVel[id]);
        }
    }

}

