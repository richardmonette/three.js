/**
 * @author mrdoob / http://mrdoob.com/
 */

import { Color } from '../../math/Color.js';
import { Matrix4 } from '../../math/Matrix4.js';
import { Vector2 } from '../../math/Vector2.js';
import { Vector3 } from '../../math/Vector3.js';

function UniformsCache() {

	var lights = {};

	return {

		get: function ( light ) {

			if ( lights[ light.id ] !== undefined ) {

				return lights[ light.id ];

			}

			var uniforms;

			switch ( light.type ) {

				case 'ProbeLight':
					uniforms = {
						position: new Vector3(),
						coeff0: new Color(),
						coeff1: new Color(),
						coeff2: new Color(),
						coeff3: new Color(),
						coeff4: new Color(),
						coeff5: new Color(),
						coeff6: new Color(),
						coeff7: new Color(),
						coeff8: new Color()
					};
					break;

				case 'DirectionalLight':
					uniforms = {
						direction: new Vector3(),
						color: new Color(),

						shadow: false,
						shadowBias: 0,
						shadowRadius: 1,
						shadowMapSize: new Vector2()
					};
					break;

				case 'SpotLight':
					uniforms = {
						position: new Vector3(),
						direction: new Vector3(),
						color: new Color(),
						distance: 0,
						coneCos: 0,
						penumbraCos: 0,
						decay: 0,

						shadow: false,
						shadowBias: 0,
						shadowRadius: 1,
						shadowMapSize: new Vector2()
					};
					break;

				case 'PointLight':
					uniforms = {
						position: new Vector3(),
						color: new Color(),
						distance: 0,
						decay: 0,

						shadow: false,
						shadowBias: 0,
						shadowRadius: 1,
						shadowMapSize: new Vector2(),
						shadowCameraNear: 1,
						shadowCameraFar: 1000
					};
					break;

				case 'HemisphereLight':
					uniforms = {
						direction: new Vector3(),
						skyColor: new Color(),
						groundColor: new Color()
					};
					break;

				case 'RectAreaLight':
					uniforms = {
						color: new Color(),
						position: new Vector3(),
						halfWidth: new Vector3(),
						halfHeight: new Vector3()
						// TODO (abelnation): set RectAreaLight shadow uniforms
					};
					break;

			}

			lights[ light.id ] = uniforms;

			return uniforms;

		}

	};

}

var count = 0;

function WebGLLights() {

	var cache = new UniformsCache();

	var state = {

		id: count ++,

		hash: {
			stateID: - 1,
			probeLength: -1,
			directionalLength: - 1,
			pointLength: - 1,
			spotLength: - 1,
			rectAreaLength: - 1,
			hemiLength: - 1,
			shadowsLength: - 1
		},

		ambient: [ 0, 0, 0 ],
		probe: [],
		directional: [],
		directionalShadowMap: [],
		directionalShadowMatrix: [],
		spot: [],
		spotShadowMap: [],
		spotShadowMatrix: [],
		rectArea: [],
		point: [],
		pointShadowMap: [],
		pointShadowMatrix: [],
		hemi: []

	};

	var vector3 = new Vector3();
	var matrix4 = new Matrix4();
	var matrix42 = new Matrix4();

	function setup( lights, shadows, camera ) {

		var r = 0, g = 0, b = 0;

		var probeLength = 0;
		var directionalLength = 0;
		var pointLength = 0;
		var spotLength = 0;
		var rectAreaLength = 0;
		var hemiLength = 0;

		var viewMatrix = camera.matrixWorldInverse;

		for ( var i = 0, l = lights.length; i < l; i ++ ) {

			var light = lights[ i ];

			var color = light.color;
			var intensity = light.intensity;
			var distance = light.distance;

			var shadowMap = ( light.shadow && light.shadow.map ) ? light.shadow.map.texture : null;

			if ( light.isAmbientLight ) {

				r += color.r * intensity;
				g += color.g * intensity;
				b += color.b * intensity;

			} else if ( light.isProbeLight ) {

				var uniforms = cache.get( light );

				uniforms.position.setFromMatrixPosition( light.matrixWorld );
				uniforms.position.applyMatrix4( viewMatrix );

				uniforms.coeff0.copy( light.coefficients[0] ).multiplyScalar( intensity );
				uniforms.coeff1.copy( light.coefficients[1] ).multiplyScalar( intensity );
				uniforms.coeff2.copy( light.coefficients[2] ).multiplyScalar( intensity );
				uniforms.coeff3.copy( light.coefficients[3] ).multiplyScalar( intensity );
				uniforms.coeff4.copy( light.coefficients[4] ).multiplyScalar( intensity );
				uniforms.coeff5.copy( light.coefficients[5] ).multiplyScalar( intensity );
				uniforms.coeff6.copy( light.coefficients[6] ).multiplyScalar( intensity );
				uniforms.coeff7.copy( light.coefficients[7] ).multiplyScalar( intensity );
				uniforms.coeff8.copy( light.coefficients[8] ).multiplyScalar( intensity );

				state.probe[ probeLength ] = uniforms;

				probeLength ++;

			} else if ( light.isDirectionalLight ) {

				var uniforms = cache.get( light );

				uniforms.color.copy( light.color ).multiplyScalar( light.intensity );
				uniforms.direction.setFromMatrixPosition( light.matrixWorld );
				vector3.setFromMatrixPosition( light.target.matrixWorld );
				uniforms.direction.sub( vector3 );
				uniforms.direction.transformDirection( viewMatrix );

				uniforms.shadow = light.castShadow;

				if ( light.castShadow ) {

					var shadow = light.shadow;

					uniforms.shadowBias = shadow.bias;
					uniforms.shadowRadius = shadow.radius;
					uniforms.shadowMapSize = shadow.mapSize;

				}

				state.directionalShadowMap[ directionalLength ] = shadowMap;
				state.directionalShadowMatrix[ directionalLength ] = light.shadow.matrix;
				state.directional[ directionalLength ] = uniforms;

				directionalLength ++;

			} else if ( light.isSpotLight ) {

				var uniforms = cache.get( light );

				uniforms.position.setFromMatrixPosition( light.matrixWorld );
				uniforms.position.applyMatrix4( viewMatrix );

				uniforms.color.copy( color ).multiplyScalar( intensity );
				uniforms.distance = distance;

				uniforms.direction.setFromMatrixPosition( light.matrixWorld );
				vector3.setFromMatrixPosition( light.target.matrixWorld );
				uniforms.direction.sub( vector3 );
				uniforms.direction.transformDirection( viewMatrix );

				uniforms.coneCos = Math.cos( light.angle );
				uniforms.penumbraCos = Math.cos( light.angle * ( 1 - light.penumbra ) );
				uniforms.decay = light.decay;

				uniforms.shadow = light.castShadow;

				if ( light.castShadow ) {

					var shadow = light.shadow;

					uniforms.shadowBias = shadow.bias;
					uniforms.shadowRadius = shadow.radius;
					uniforms.shadowMapSize = shadow.mapSize;

				}

				state.spotShadowMap[ spotLength ] = shadowMap;
				state.spotShadowMatrix[ spotLength ] = light.shadow.matrix;
				state.spot[ spotLength ] = uniforms;

				spotLength ++;

			} else if ( light.isRectAreaLight ) {

				var uniforms = cache.get( light );

				// (a) intensity is the total visible light emitted
				//uniforms.color.copy( color ).multiplyScalar( intensity / ( light.width * light.height * Math.PI ) );

				// (b) intensity is the brightness of the light
				uniforms.color.copy( color ).multiplyScalar( intensity );

				uniforms.position.setFromMatrixPosition( light.matrixWorld );
				uniforms.position.applyMatrix4( viewMatrix );

				// extract local rotation of light to derive width/height half vectors
				matrix42.identity();
				matrix4.copy( light.matrixWorld );
				matrix4.premultiply( viewMatrix );
				matrix42.extractRotation( matrix4 );

				uniforms.halfWidth.set( light.width * 0.5, 0.0, 0.0 );
				uniforms.halfHeight.set( 0.0, light.height * 0.5, 0.0 );

				uniforms.halfWidth.applyMatrix4( matrix42 );
				uniforms.halfHeight.applyMatrix4( matrix42 );

				// TODO (abelnation): RectAreaLight distance?
				// uniforms.distance = distance;

				state.rectArea[ rectAreaLength ] = uniforms;

				rectAreaLength ++;

			} else if ( light.isPointLight ) {

				var uniforms = cache.get( light );

				uniforms.position.setFromMatrixPosition( light.matrixWorld );
				uniforms.position.applyMatrix4( viewMatrix );

				uniforms.color.copy( light.color ).multiplyScalar( light.intensity );
				uniforms.distance = light.distance;
				uniforms.decay = light.decay;

				uniforms.shadow = light.castShadow;

				if ( light.castShadow ) {

					var shadow = light.shadow;

					uniforms.shadowBias = shadow.bias;
					uniforms.shadowRadius = shadow.radius;
					uniforms.shadowMapSize = shadow.mapSize;
					uniforms.shadowCameraNear = shadow.camera.near;
					uniforms.shadowCameraFar = shadow.camera.far;

				}

				state.pointShadowMap[ pointLength ] = shadowMap;
				state.pointShadowMatrix[ pointLength ] = light.shadow.matrix;
				state.point[ pointLength ] = uniforms;

				pointLength ++;

			} else if ( light.isHemisphereLight ) {

				var uniforms = cache.get( light );

				uniforms.direction.setFromMatrixPosition( light.matrixWorld );
				uniforms.direction.transformDirection( viewMatrix );
				uniforms.direction.normalize();

				uniforms.skyColor.copy( light.color ).multiplyScalar( intensity );
				uniforms.groundColor.copy( light.groundColor ).multiplyScalar( intensity );

				state.hemi[ hemiLength ] = uniforms;

				hemiLength ++;

			}

		}

		state.ambient[ 0 ] = r;
		state.ambient[ 1 ] = g;
		state.ambient[ 2 ] = b;

		state.probe.length = probeLength;
		state.directional.length = directionalLength;
		state.spot.length = spotLength;
		state.rectArea.length = rectAreaLength;
		state.point.length = pointLength;
		state.hemi.length = hemiLength;

		state.hash.stateID = state.id;
		state.hash.probeLength = probeLength;
		state.hash.directionalLength = directionalLength;
		state.hash.pointLength = pointLength;
		state.hash.spotLength = spotLength;
		state.hash.rectAreaLength = rectAreaLength;
		state.hash.hemiLength = hemiLength;
		state.hash.shadowsLength = shadows.length;

	}

	return {
		setup: setup,
		state: state
	};

}


export { WebGLLights };
