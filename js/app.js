import * as THREE from './three.module.js'
import { GUI } from './dat.gui.module.js'
import Engine from './engine.js'

//window.addEventListener('load', init, false);


// Constants
const boxWidth = 1.0;
const boxDepth = 0.8;
const boxHeight = 0.7;
// GUI setup
const gui = new GUI();
let initial_particles_n = 2000;
let initial_particle_radius = 0.007;
let initial_mass = 1.0;	    
let k = 120; //gas constant		
let rho0 = 0; //rest density
let mu = 5; //viscosity		
let gx = 0;	//gravity		
let gy = -9;			
let gz = 0;
let h = 1.01; //smoothing length
let initial_timestep= 1.0/60.0;
//particles
let particleMeshes = [];
const panelOptions = {
    timestep: initial_timestep,
    particleCount: initial_particles_n,
    particleRadius: initial_particle_radius,
    mass: initial_mass,
    restDensity: rho0,
    viscosity: mu,
    gasConstant: k,
    smoothingLength: h,
    gravityX: gx,
    gravityY: gy,
    gravityZ: gz,
    lightX: 1, 
    lightY: 2, 
    lightZ: 2,
    mouseForce: 1, 
    startSim: function() {
        startSimulation();
    },
    resetSim: function() {
        resetSimulation();
    }
};
//mouse events handling
let isRotating = false;
let previousMousePosition = { x: 0, y: 0 };


const halfWidth = boxWidth/2;
const halfHeight = boxHeight/2;
const halfDepth = boxDepth/2;

let engine = new Engine(boxWidth, boxHeight, boxDepth, -halfWidth, halfWidth, -halfHeight, halfHeight, -halfDepth, halfDepth);

//scene constants
const renderer = new THREE.WebGLRenderer();
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const geometry = new THREE.BoxGeometry(boxWidth, boxHeight, boxDepth);
const cubeMaterial = new THREE.MeshPhongMaterial({
    color: 0x000000,
    transparent: true,
    wireframe: true,
    wireframeLinewidth: 1,
});

const particleMaterial = new THREE.MeshPhongMaterial({ color: 0x1E90FF, transparent: true });
const cube = new THREE.Mesh(geometry, cubeMaterial);

// Lights
const ambientLight = new THREE.AmbientLight(0x404040, 2.0); // soft white light
scene.add(ambientLight);

const pointLight = new THREE.PointLight(0xffffff, 1, 100);
pointLight.position.set(panelOptions.lightX, panelOptions.lightY, panelOptions.lightZ);
scene.add(pointLight);

// Light helper sphere
const lightHelperGeometry = new THREE.SphereGeometry(0.05, 20, 10);
const lightHelperMaterial = new THREE.MeshPhongMaterial({ color: 0x000000, transparent: true, wireframe: true });
const lightHelper = new THREE.Mesh(lightHelperGeometry, lightHelperMaterial);
lightHelper.position.copy(pointLight.position);
scene.add(lightHelper);



function updateParticleCount(n) {
    particleMeshes.forEach(mesh => scene.remove(mesh));
    // Reset the engine if the new count is less than the current count
    if (n < particleMeshes.length) {
        engine.reset();
    }
    particleMeshes = new Array(n);

    let particleGeometry = new THREE.SphereGeometry(panelOptions.particleRadius, 16, 8);
    for (var i = 0; i < n; i++) {
        var m = new THREE.Mesh(particleGeometry, particleMaterial);
        m.position.x = (Math.random() - 0.5) * boxWidth;
        m.position.y = (Math.random() * 0.5) * boxHeight; //put particles from half to top
        m.position.z = (Math.random() - 0.5) * boxDepth;
        //console.log(" [SET NUM PARTICLES APP] Particle MESH after definition ", m);
        particleMeshes[i] = m;
        //console.log("[SET NUM PARTICLES APP] Particle meshes after setting i-th array entry ", i, particleMeshes);
        scene.add(m);
    }
    //console.log("[SET NUM PARTICLES APP] Particle meshes after definition ", particleMeshes);
    engine.updateParticleCount(n, particleMeshes, panelOptions.mass);
    renderer.render(scene, camera);
}


function init(){
    updateSimulationParameters();
    initScene();
    addListeners();
    updateParticleCount(panelOptions.particleCount);
    setGuiPanel();
}


function initScene(){
    camera.position.z = 1;
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setClearColor(0xFAEBD7, 1);  // Set background color to antiquewhite
    document.body.appendChild(renderer.domElement); //renderer creates a canvas element that is used to draw and render the scene
    scene.add(cube);
}

function updateSimulationParameters() {
    let h= panelOptions.smoothingLength;
    let mass = panelOptions.mass;
    let gasConstant = panelOptions.gasConstant;
    let restDensity = panelOptions.restDensity;
    let viscosity = panelOptions.viscosity;
    let ts= panelOptions.timestep;
    let mouseMultiplier= panelOptions.mouseForce;
    let gr_x = panelOptions.gravityX;
    let gr_y = panelOptions.gravityY;
    let gr_z = panelOptions.gravityZ;
    engine.setSimulationParameters(mass, gasConstant, restDensity, viscosity, h, ts, mouseMultiplier, gr_x, gr_y, gr_z);
}

function updateLightPosition() {
    pointLight.position.set(panelOptions.lightX, panelOptions.lightY, panelOptions.lightZ);
    lightHelper.position.copy(pointLight.position);
}


function setGuiPanel(){
    gui.add(panelOptions, 'particleCount', 0, 10000).name('Particle count').step(1).onChange(value => {
        updateParticleCount(panelOptions.particleCount);
    });
    
    gui.add(panelOptions, 'particleRadius', 0.001, 0.020).name('Particle radius').step(0.001).onChange(value => {
        panelOptions.particleRadius = value;
        updateParticleCount(panelOptions.particleCount);
    });
    gui.add(panelOptions, 'mass', 1, 3.0).name('Particle mass').step(0.01).onChange(updateSimulationParameters);
    gui.add(panelOptions, 'gasConstant', 1, 1000).name('Gas constant').onChange(updateSimulationParameters);
    gui.add(panelOptions, 'restDensity', 0, 5).step(0.1).name('Rest density').onChange(updateSimulationParameters);
    gui.add(panelOptions, 'viscosity', 0, 11).name('Viscosity').onChange(updateSimulationParameters);
    gui.add(panelOptions, 'smoothingLength', 1, 1.25).name('Smoothing length').step(0.001).onChange(updateSimulationParameters);
    gui.add(panelOptions, 'gravityX', -100, 100).step(1).name('Gravity X').onChange(updateSimulationParameters);
    gui.add(panelOptions, 'gravityY', -100, 100).step(1).name('Gravity Y').onChange(updateSimulationParameters);
    gui.add(panelOptions, 'gravityZ', -100, 100).step(1).name('Gravity Z').onChange(updateSimulationParameters);
    gui.add(panelOptions, 'lightX', -10, 10).step(0.01).name('Light X').onChange(updateLightPosition);
    gui.add(panelOptions, 'lightY', -10, 10).step(0.01).name('Light Y').onChange(updateLightPosition);
    gui.add(panelOptions, 'lightZ', -10, 10).step(0.01).name('Light Z').onChange(updateLightPosition);
    gui.add(panelOptions, 'mouseForce', 0, 5).name('Mouse multiplier').step(0.01).onChange(updateSimulationParameters);
    gui.add(panelOptions, 'timestep', 0.01, 0.021).name('Time step').step(0.001).onChange(updateSimulationParameters);

    gui.add(panelOptions, 'startSim').name('Start Simulation');
    gui.add(panelOptions, 'resetSim').name('Reset Simulation');
}

function getMousePosition(event, element){
    var boundingRect = element.getBoundingClientRect();
    return {
        //event.clientX provides the horizontal coordinate of the mouse cursor relative to the viewport.
        //boundingRect.left is the distance from the left edge of the viewport to the left edge of the element.
        x: event.clientX - boundingRect.left,
        y: event.clientY - boundingRect.top
    };
}

function mouseDownHandler(event) {
    isRotating = true;
    previousMousePosition = { x: event.clientX, y: event.clientY };
}

//to force velocity
function castRay(x, y, canvas, camera) {
    // get the normalized device coordinates (NDC) of the mouse, required for setFromCamera API of raycaster
    let ndcX = (x / canvas.clientWidth) * 2 - 1;
    let ndcY = - (y / canvas.clientHeight) * 2 + 1;  // "-" sign is due to the fact that canvas y-coordinate increases downward while 
    //the NDC y-coordinate increases upward, so the conversion needs to invert the direction
    
    // Create a ray from the camera through the mouse position in NDC
    let raycaster = new THREE.Raycaster();
    let mouseVector = new THREE.Vector2(ndcX, ndcY);
    raycaster.setFromCamera(mouseVector, camera);

    const intersects = raycaster.intersectObjects( scene.children );
    if(intersects.length == 0) return null;
    else return intersects;
}

function mouseMoveHandler(event) {
    if (isRotating) {
        let deltaMove = { x: event.clientX - previousMousePosition.x, y: event.clientY - previousMousePosition.y };
        console.log(deltaMove);
        // rotate camera around the cube
        const rotationSpeed = 0.005;
        const spherical = new THREE.Spherical();  // a point spherical coordinates
        spherical.setFromVector3(camera.position);
        //adjust spherical coordinates, where theta is the horizontal angle and phi is the vertical angle
        spherical.theta -= deltaMove.x * rotationSpeed;  //azimuthal angle  (horizontal rotation)
        spherical.phi -= deltaMove.y * rotationSpeed;  //polar angle (vertical rotation)
        //reconvert into cartesian coordinates
        camera.position.setFromSpherical(spherical);
        //ensure camera always looks at the center of the cube
        camera.lookAt(cube.position);
        previousMousePosition = { x: event.clientX, y: event.clientY };
        renderer.render(scene, camera);
    }

    let mousePos = getMousePosition(event, renderer.domElement);
    let int_Objects = castRay(mousePos.x, mousePos.y, renderer.domElement, camera);
    //console.log("[FORCE VELOCITY] intersected objects ", int_Objects);
    if(int_Objects != null) {
        engine.applyMouseForce(int_Objects, event.movementX, event.movementY);
    }
}


function resetSimulation() {
    // Reset the engine and scene
    engine = new Engine(boxWidth, boxHeight, boxDepth, -halfWidth, halfWidth, -halfHeight, halfHeight, -halfDepth, halfDepth);

    particleMeshes.forEach(mesh => scene.remove(mesh));
    particleMeshes = [];
    
    // Reinitialize the scene
    initScene();
    updateParticleCount(panelOptions.particleCount);
    updateSimulationParameters();
    updateLightPosition();
    renderer.render(scene, camera);
}


function mouseWheelHandler(event) {
    let zoomFactor = 1.05; // Adjust this factor to control the zoom speed
    if (event.deltaY < 0) {
       camera.fov /= zoomFactor; // Zoom in
    } else {
        camera.fov *= zoomFactor; // Zoom out
    }
    camera.updateProjectionMatrix();
    renderer.render(scene, camera);
}
    
  

function mouseUpHandler(event) {
    isRotating = false;
}

function addListeners(){
    //mouse event listeners
    renderer.domElement.addEventListener('mousedown', mouseDownHandler, false);
    renderer.domElement.addEventListener('mousemove', mouseMoveHandler, false);
    renderer.domElement.addEventListener('mouseup', mouseUpHandler, false);
    renderer.domElement.addEventListener('wheel', mouseWheelHandler, false);
}

let clock = 0;

function startSimulation() {
    requestAnimationFrame(startSimulation);
    //if(clock < 7){
        engine.simulationStep();
        for (var i = 0; i < particleMeshes.length; i++){
            //console.log("[ANIMATE] particleMeshes before update is ", particleMeshes);
            engine.updateMeshPosition(i, particleMeshes[i].position);
            //console.log("[ANIMATE] particleMeshes after update is ", particleMeshes);
        }
        renderer.render(scene, camera);
    //}
    //clock++

}
init();